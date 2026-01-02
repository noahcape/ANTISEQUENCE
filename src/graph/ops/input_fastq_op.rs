use needletail::*;

use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::Arc;
use std::io::Write;
use parking_lot::Mutex;

use crate::errors::*;
use crate::expr::LabelOrAttr;
use crate::graph::*;
use std::sync::OnceLock;

/// Returns the chunk size for batch processing, lazily initialized.
///
/// The chunk size determines how many reads are processed together in a single batch.
/// This affects memory usage and cache locality.
///
/// It checks the `ANTISEQ_CHUNK_SIZE` environment variable.
/// - If set to a valid positive integer, that value is used.
/// - Defaults to 512 if unset or invalid.
///
/// # Example
/// ```sh
/// export ANTISEQ_CHUNK_SIZE=1024
/// ```
fn chunk_size() -> usize {
    static CHUNK: OnceLock<usize> = OnceLock::new();
    *CHUNK.get_or_init(|| {
        std::env::var("ANTISEQ_CHUNK_SIZE")
            .ok()
            .and_then(|s| s.parse::<usize>().ok())
            .filter(|&v| v > 0)
            .unwrap_or(512)
    })
}

/// Progress reporting interval (number of reads between updates).
/// Defaults to 1,000,000 reads. Set SEQPROC_PROGRESS_INTERVAL to customize.
fn progress_interval() -> usize {
    static INTERVAL: OnceLock<usize> = OnceLock::new();
    *INTERVAL.get_or_init(|| {
        std::env::var("SEQPROC_PROGRESS_INTERVAL")
            .ok()
            .and_then(|s| s.parse::<usize>().ok())
            .filter(|&v| v > 0)
            .unwrap_or(1_000_000)
    })
}

/// Check if progress reporting is enabled. Set SEQPROC_QUIET=1 to disable.
fn progress_enabled() -> bool {
    static ENABLED: OnceLock<bool> = OnceLock::new();
    *ENABLED.get_or_init(|| {
        !std::env::var("SEQPROC_QUIET")
            .ok()
            .map(|v| v == "1" || v.eq_ignore_ascii_case("true"))
            .unwrap_or(false)
    })
}

/// An input operation that reads FASTQ records and emits batches of `Read`s.
///
/// This struct manages shared access to one or more FASTQ readers (e.g., from files or stdin).
/// It is designed to be thread-safe, allowing multiple worker threads to pull batches of reads
/// concurrently from the same source(s).
///
/// # Key Features
/// - Supports reading from single or multiple files/readers.
/// - Supports interleaved paired-end reads (where R1 and R2 are in the same file).
/// - Tracks statistics like read counts and length distributions.
pub struct InputFastqOp<'reader> {
    readers: Vec<(Mutex<Box<dyn FastxReader + 'reader>>, Arc<Origin>)>,
    idx: AtomicUsize,
    interleaved: usize,
    read_counts: Vec<AtomicUsize>,
    read_length_min: Vec<AtomicUsize>,
    read_length_max: Vec<AtomicUsize>,
    read_length_sum: Vec<AtomicUsize>,
    /// Tracks the last progress milestone reported (in units of progress_interval)
    last_progress_milestone: AtomicUsize,
}

impl<'reader> InputFastqOp<'reader> {
    const NAME: &'static str = "InputFastqOp";

    /// Stream reads created from fastq records from an input file.
    pub fn from_file(file: impl AsRef<str>) -> Result<Self> {
        let reader = Mutex::new(parse_fastx_file(file.as_ref()).map_err(|e| Error::FileIo {
            file: file.as_ref().to_owned(),
            source: Box::new(e),
        })?);
        let n_fastqs = 1;       // Only track one fastq stream by default

        // Initialize the actual InputFastqOp struct since we prepped the file Mutex and started the 
        // statistics tracking vectors
        Ok(Self {
            readers: vec![(reader, Arc::new(Origin::File(file.as_ref().to_owned())))],          // Number via atomic ref counting for thread safety
            idx: AtomicUsize::new(0),                                                           // Index of the current reader
            interleaved: 1,                                                                     // Interleaved paired-end reads
            read_counts: (0..n_fastqs).map(|_| AtomicUsize::new(0)).collect(),                  // Read counts
            read_length_min: (0..n_fastqs).map(|_| AtomicUsize::new(usize::MAX)).collect(),     // Minimum read length clamped to usize::MAX to ensure we find the minimum
            read_length_max: (0..n_fastqs).map(|_| AtomicUsize::new(0)).collect(),              // Maximum read length clamped to 0 to ensure we find the maximum
            read_length_sum: (0..n_fastqs).map(|_| AtomicUsize::new(0)).collect(),              // Sum of read lengths
            last_progress_milestone: AtomicUsize::new(0),
        })
    }

    /// Stream reads created from fastq records from multiple input files.
    pub fn from_files<S: AsRef<str>>(files: impl IntoIterator<Item = S>) -> Result<Self> {
        let readers = files
            .into_iter()
            .map(|f| {
                let file = f.as_ref();
                (
                    Mutex::new(parse_fastx_file(file).unwrap_or_else(|e| panic!("{e}"))),
                    Arc::new(Origin::File(file.to_owned())),
                )
            })
            .collect::<Vec<_>>();

        let n_fastqs = readers.len();

        Ok(Self {
            readers,
            idx: AtomicUsize::new(0),
            interleaved: 1,
            read_counts: (0..n_fastqs).map(|_| AtomicUsize::new(0)).collect(),
            read_length_min: (0..n_fastqs).map(|_| AtomicUsize::new(usize::MAX)).collect(),
            read_length_max: (0..n_fastqs).map(|_| AtomicUsize::new(0)).collect(),
            read_length_sum: (0..n_fastqs).map(|_| AtomicUsize::new(0)).collect(),
            last_progress_milestone: AtomicUsize::new(0),
        })
    }

    /// Stream reads created from interleaved fastq records from an input file.
    pub fn from_file_interleaved(file: impl AsRef<str>, interleaved: usize) -> Result<Self> {
        let reader = Mutex::new(parse_fastx_file(file.as_ref()).map_err(|e| Error::FileIo {
            file: file.as_ref().to_owned(),
            source: Box::new(e),
        })?);
        let n_fastqs = interleaved;

        Ok(Self {
            readers: vec![(reader, Arc::new(Origin::File(file.as_ref().to_owned())))],
            idx: AtomicUsize::new(0),
            interleaved,
            read_counts: (0..n_fastqs).map(|_| AtomicUsize::new(0)).collect(),
            read_length_min: (0..n_fastqs).map(|_| AtomicUsize::new(usize::MAX)).collect(),
            read_length_max: (0..n_fastqs).map(|_| AtomicUsize::new(0)).collect(),
            read_length_sum: (0..n_fastqs).map(|_| AtomicUsize::new(0)).collect(),
            last_progress_milestone: AtomicUsize::new(0),
        })
    }

    /// Stream reads created from fastq records from an arbitrary `Read`er.
    pub fn from_reader(reader: impl std::io::Read + Send + 'reader) -> Result<Self> {
        let reader =
            Mutex::new(parse_fastx_reader(reader).map_err(|e| Error::BytesIo(Box::new(e)))?);
        let n_fastqs = 1;

        Ok(Self {
            readers: vec![(reader, Arc::new(Origin::Bytes))],
            idx: AtomicUsize::new(0),
            interleaved: 1,
            read_counts: (0..n_fastqs).map(|_| AtomicUsize::new(0)).collect(),
            read_length_min: (0..n_fastqs).map(|_| AtomicUsize::new(usize::MAX)).collect(),
            read_length_max: (0..n_fastqs).map(|_| AtomicUsize::new(0)).collect(),
            read_length_sum: (0..n_fastqs).map(|_| AtomicUsize::new(0)).collect(),
            last_progress_milestone: AtomicUsize::new(0),
        })
    }

    /// Stream reads created from fastq records from multiple arbitrary `Read`ers.
    pub fn from_readers<R: std::io::Read + Send + 'reader>(
        readers: impl IntoIterator<Item = R>,
    ) -> Result<Self> {
        let readers = readers
            .into_iter()
            .map(|r| {
                (
                    Mutex::new(parse_fastx_reader(r).unwrap_or_else(|e| panic!("{e}"))),
                    Arc::new(Origin::Bytes),
                )
            })
            .collect::<Vec<_>>();

        let n_fastqs = readers.len();

        Ok(Self {
            readers,
            idx: AtomicUsize::new(0),
            interleaved: 1,
            read_counts: (0..n_fastqs).map(|_| AtomicUsize::new(0)).collect(),
            read_length_min: (0..n_fastqs).map(|_| AtomicUsize::new(usize::MAX)).collect(),
            read_length_max: (0..n_fastqs).map(|_| AtomicUsize::new(0)).collect(),
            read_length_sum: (0..n_fastqs).map(|_| AtomicUsize::new(0)).collect(),
            last_progress_milestone: AtomicUsize::new(0),
        })
    }

    /// Stream reads created from interleaved fastq records from an arbitrary `Read`er.
    pub fn from_interleaved_reader(
        reader: impl std::io::Read + Send + 'reader,
        interleaved: usize,
    ) -> Result<Self> {
        let reader =
            Mutex::new(parse_fastx_reader(reader).map_err(|e| Error::BytesIo(Box::new(e)))?);
        let n_fastqs = interleaved;

        Ok(Self {
            readers: vec![(reader, Arc::new(Origin::Bytes))],
            idx: AtomicUsize::new(0),
            interleaved,
            read_counts: (0..n_fastqs).map(|_| AtomicUsize::new(0)).collect(),
            read_length_min: (0..n_fastqs).map(|_| AtomicUsize::new(usize::MAX)).collect(),
            read_length_max: (0..n_fastqs).map(|_| AtomicUsize::new(0)).collect(),
            read_length_sum: (0..n_fastqs).map(|_| AtomicUsize::new(0)).collect(),
            last_progress_milestone: AtomicUsize::new(0),
        })
    }
}

impl<'reader, T: Trace> GraphNode<T> for InputFastqOp<'reader> {
    /// Main execution logic for the InputOp.
    ///
    /// This method is called by the graph executor. It produces a batch of reads by:
    /// 1. locking the readers (amortized over the chunk size).
    /// 2. parsing FASTQ records.
    /// 3. populating `Read` objects, reusing memory where possible.
    ///
    /// # Arguments
    /// * `reads` - An optional vector of `Read`s to reuse (recycle). If `None`, a new vector is allocated.
    /// * `trace` - Tracing infrastructure for performance monitoring.
    ///
    /// # Returns
    /// * `Ok((Some(batch), false))` - A batch of reads was successfully read.
    /// * `Ok((None, true))` - End of input reached (EOF).
    fn run(&self, reads: Option<Vec<Read>>, trace: &T) -> Result<(Option<Vec<Read>>, bool)> {
        let start = trace.start(&reads);
        let cs = chunk_size();
        // If a vector of reads was passed in (from a downstream op recycling it), use it.
        // Otherwise, allocate a new vector with capacity for the chunk size.
        // We reduce a ton of churn here!
        let mut b = reads.unwrap_or_else(|| Vec::with_capacity(cs));
        // Do NOT clear b here, we want to reuse its elements (Read objects and their buffers).
        // This gets us a massive speedup :)

        // Lock readers once per batch refill to amortize lock acquisition overhead.
        // This is critical for performance when multiple threads contend for the same readers.
        let mut locked_readers = self
            .readers
            .iter()
            // r is a MutexGuard, o is an Arc<Origin> holding the origin of the reader
            .map(|(r, o)| (r.lock(), o))
            .collect::<Vec<_>>();

        let mut i = 0;
        // outer loop over the chunk size
        'outer: for _ in 0..cs {
            // Atomically increment global index to assign unique IDs to reads across threads.
            let idx = self.idx.fetch_add(self.interleaved, Ordering::Relaxed);
            
            if self.interleaved > 1 {
                // CASE: Interleaved input (e.g., R1 and R2 in same file).
                // All records come from the first reader.
                let (locked_reader, origin) = &mut locked_readers[0];

                // If the batch vector needs more Read objects, create new ones.
                if i >= b.len() {
                    b.push(Read::new());
                }
                let curr_read = &mut b[i];
                // Note: We do NOT call curr_read.clear(). `set_fastq_entry` handles overwriting/recycling.

                let mut slot_idx = 0;
                
                // j here is the index of the specific read within the current group.
                // For standard paired-end sequencing this is 2.
                // Thus in that 2 case, j=0 represents read 1 (forward read)
                // and j=1 represents read 2 (reverse read).
                // So this loop bundles the two reads into a single Read object.
                for j in 0..self.interleaved {
                    // Fetch the next record from the locked reader.
                    let Some(record) = locked_reader.next() else {
                        if j == 0 {
                            // EOF reached cleanly between pairs.
                            b.truncate(i); // Remove any unused stale reads at the end of the batch.
                            break 'outer;
                        }
                        // EOF reached in the middle of a pair/tuple -> Error.
                        Err(Error::UnpairedRead(format!("\"{}\"", &**origin)))?
                    };
                    let record = record.map_err(|e| Error::ParseRecord {
                        origin: (***origin).clone(),
                        idx: idx + j,
                        source: Box::new(e),
                    })?;
                    
                    let seq = record.seq();
                    self.update_length_stats(j, seq.len());
                    
                    // Populate Name entry (e.g., @read1/1)
                    curr_read.set_fastq_entry(
                        slot_idx,
                        StrType::Name((j + 1) as _),
                        record.id(),
                        None,
                        Arc::clone(origin),
                        idx + j,
                    );
                    slot_idx += 1;
                    
                    // Populate Sequence entry (and Qual)
                    curr_read.set_fastq_entry(
                        slot_idx,
                        StrType::Seq((j + 1) as _),
                        &seq,
                        record.qual(),
                        Arc::clone(origin),
                        idx + j,
                    );
                    slot_idx += 1;
                }
            } else {
                // CASE: Separate files (e.g., R1 file and R2 file).
                // We iterate over all readers in lock-step.
                if i >= b.len() {
                    b.push(Read::new());
                }
                let curr_read = &mut b[i];

                // Loop through every open file reader. j is the index of the file (0 for R1, 1 for R2, etc.)
                // For every iteration of the outer loop (which represents one "read tuple"), 
                // this inner loop pulls one record from every file.
                let mut slot_idx = 0;
                for (j, (locked_reader, origin)) in locked_readers.iter_mut().enumerate() {
                    let Some(record) = locked_reader.next() else {
                        if j == 0 {
                            b.truncate(i);
                            break 'outer;           // EOF reached cleanly between pairs.
                        }
                        Err(Error::UnpairedRead(format!("\"{}\"", &**origin)))?
                    };
                    let record = record.map_err(|e| Error::ParseRecord {
                        origin: (***origin).clone(),
                        idx,
                        source: Box::new(e),
                    })?;
                    let seq = record.seq();
                    self.update_length_stats(j, seq.len());
                    
                    // combine the separate records into a single Read (cur_read) object
                    curr_read.set_fastq_entry(
                        slot_idx,
                        StrType::Name((j + 1) as _),
                        record.id(),
                        None,
                        Arc::clone(origin),
                        idx,
                    );
                    slot_idx += 1;
                    
                    curr_read.set_fastq_entry(
                        slot_idx,
                        StrType::Seq((j + 1) as _),
                        &seq,
                        record.qual(),
                        Arc::clone(origin),
                        idx,
                    );
                    slot_idx += 1;
                }
            }
            i += 1;
        }
        
        // Trim vector to actual number of reads read (handles EOF partial batch).
        if b.len() > i {
            b.truncate(i);
        }

        // If batch is empty, we hit EOF immediately. Were done!
        if b.is_empty() {
            return Ok((None, true));
        }

        // Progress reporting: check if we crossed a milestone
        if progress_enabled() {
            let interval = progress_interval();
            let current_idx = self.idx.load(Ordering::Relaxed);
            let current_milestone = current_idx / interval;
            let last_milestone = self.last_progress_milestone.load(Ordering::Relaxed);
            
            if current_milestone > last_milestone {
                // Try to claim this milestone (avoid duplicate prints from multiple threads)
                if self.last_progress_milestone.compare_exchange(
                    last_milestone,
                    current_milestone,
                    Ordering::Relaxed,
                    Ordering::Relaxed,
                ).is_ok() {
                    let _ = write!(std::io::stderr(), "\r[seqproc] Processed {} reads...", current_milestone * interval);
                    let _ = std::io::stderr().flush();
                }
            }
        }

        let res_opt = Some(b);
        trace.add(<Self as GraphNode<T>>::name(self), start, &res_opt);
        Ok((res_opt, false))
    }

    fn required_names(&self) -> &[LabelOrAttr] {
        &[]
    }

    fn name(&self) -> &'static str {
        Self::NAME
    }

    fn input_stats(&self) -> Option<InputStats> {
        let n_fastqs = self.read_counts.len();

        let mut read_counts = Vec::with_capacity(n_fastqs);
        let mut read_length_min = Vec::with_capacity(n_fastqs);
        let mut read_length_max = Vec::with_capacity(n_fastqs);
        let mut read_length_sum = Vec::with_capacity(n_fastqs);

        for i in 0..n_fastqs {
            let count = self.read_counts[i].load(Ordering::Relaxed);
            read_counts.push(count);

            let min = self.read_length_min[i].load(Ordering::Relaxed);
            let max = self.read_length_max[i].load(Ordering::Relaxed);
            let sum = self.read_length_sum[i].load(Ordering::Relaxed);

            if count == 0 {
                read_length_min.push(0);
                read_length_max.push(0);
            } else {
                read_length_min.push(min);
                read_length_max.push(max);
            }
            read_length_sum.push(sum);
        }

        Some(InputStats {
            n_fastqs,
            read_counts,
            read_length_min,
            read_length_max,
            read_length_sum,
        })
    }
}

impl<'reader> InputFastqOp<'reader> {
    fn update_length_stats(&self, lane: usize, len: usize) {
        self.read_counts[lane].fetch_add(1, Ordering::Relaxed);
        self.read_length_sum[lane].fetch_add(len, Ordering::Relaxed);

        let min = &self.read_length_min[lane];
        let mut current_min = min.load(Ordering::Relaxed);
        while len < current_min {
            match min.compare_exchange(current_min, len, Ordering::Relaxed, Ordering::Relaxed) {
                Ok(_) => break,
                Err(c) => current_min = c,
            }
        }

        let max = &self.read_length_max[lane];
        let mut current_max = max.load(Ordering::Relaxed);
        while len > current_max {
            match max.compare_exchange(current_max, len, Ordering::Relaxed, Ordering::Relaxed) {
                Ok(_) => break,
                Err(c) => current_max = c,
            }
        }
    }
}
