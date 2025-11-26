use needletail::*;

use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::Arc;
use parking_lot::Mutex;

use crate::errors::*;
use crate::expr::LabelOrAttr;
use crate::graph::*;
use std::sync::OnceLock;

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

pub struct InputFastqOp<'reader> {
    readers: Vec<(Mutex<Box<dyn FastxReader + 'reader>>, Arc<Origin>)>,
    idx: AtomicUsize,
    interleaved: usize,
    read_counts: Vec<AtomicUsize>,
    read_length_min: Vec<AtomicUsize>,
    read_length_max: Vec<AtomicUsize>,
    read_length_sum: Vec<AtomicUsize>,
}

impl<'reader> InputFastqOp<'reader> {
    const NAME: &'static str = "InputFastqOp";

    /// Stream reads created from fastq records from an input file.
    pub fn from_file(file: impl AsRef<str>) -> Result<Self> {
        let reader = Mutex::new(parse_fastx_file(file.as_ref()).map_err(|e| Error::FileIo {
            file: file.as_ref().to_owned(),
            source: Box::new(e),
        })?);
        let n_fastqs = 1;

        Ok(Self {
            readers: vec![(reader, Arc::new(Origin::File(file.as_ref().to_owned())))],
            idx: AtomicUsize::new(0),
            interleaved: 1,
            read_counts: (0..n_fastqs).map(|_| AtomicUsize::new(0)).collect(),
            read_length_min: (0..n_fastqs).map(|_| AtomicUsize::new(usize::MAX)).collect(),
            read_length_max: (0..n_fastqs).map(|_| AtomicUsize::new(0)).collect(),
            read_length_sum: (0..n_fastqs).map(|_| AtomicUsize::new(0)).collect(),
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
        })
    }
}

impl<'reader, T: Trace> GraphNode<T> for InputFastqOp<'reader> {
    fn run(&self, reads: Option<Vec<Read>>, trace: &T) -> Result<(Option<Vec<Read>>, bool)> {
        let start = trace.start(&reads);
        
        let cs = chunk_size();
        let mut b = reads.unwrap_or_else(|| Vec::with_capacity(cs));
        // Do NOT clear b here, we want to reuse its elements.

        // Lock readers once per refill to amortize lock overhead
        let mut locked_readers = self
            .readers
            .iter()
            .map(|(r, o)| (r.lock(), o))
            .collect::<Vec<_>>();

        let mut i = 0;
        'outer: for _ in 0..cs {
            let idx = self.idx.fetch_add(self.interleaved, Ordering::Relaxed);
            
            if self.interleaved > 1 {
                // interleaved records all come from one file
                let (locked_reader, origin) = &mut locked_readers[0];

                if i >= b.len() {
                    b.push(Read::new());
                }
                let curr_read = &mut b[i];
                // Removed curr_read.clear() to allow recycling

                let mut slot_idx = 0;
                for j in 0..self.interleaved {
                    let Some(record) = locked_reader.next() else {
                        if j == 0 {
                            b.truncate(i); // Remove unused reads at the end
                            break 'outer;
                        }
                        Err(Error::UnpairedRead(format!("\"{}\"", &**origin)))?
                    };
                    let record = record.map_err(|e| Error::ParseRecord {
                        origin: (***origin).clone(),
                        idx: idx + j,
                        source: Box::new(e),
                    })?;
                    
                    let seq = record.seq();
                    self.update_length_stats(j, seq.len());
                    
                    curr_read.set_fastq_entry(
                        slot_idx,
                        StrType::Name((j + 1) as _),
                        record.id(),
                        None,
                        Arc::clone(origin),
                        idx + j,
                    );
                    slot_idx += 1;
                    
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
                // gather records from multiple different files
                if i >= b.len() {
                    b.push(Read::new());
                }
                let curr_read = &mut b[i];
                // Removed curr_read.clear()

                let mut slot_idx = 0;
                for (j, (locked_reader, origin)) in locked_readers.iter_mut().enumerate() {
                    let Some(record) = locked_reader.next() else {
                        if j == 0 {
                            b.truncate(i);
                            break 'outer;
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
        
        if b.len() > i {
            b.truncate(i);
        }

        if b.is_empty() {
            return Ok((None, true));
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
