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
            .unwrap_or(256)
    })
}

pub struct InputFastqOp<'reader> {
    readers: Vec<(Mutex<Box<dyn FastxReader + 'reader>>, Arc<Origin>)>,
    idx: AtomicUsize,
    interleaved: usize,
}

impl<'reader> InputFastqOp<'reader> {
    const NAME: &'static str = "InputFastqOp";

    /// Stream reads created from fastq records from an input file.
    pub fn from_file(file: impl AsRef<str>) -> Result<Self> {
        let reader = Mutex::new(parse_fastx_file(file.as_ref()).map_err(|e| Error::FileIo {
            file: file.as_ref().to_owned(),
            source: Box::new(e),
        })?);

        Ok(Self {
            readers: vec![(reader, Arc::new(Origin::File(file.as_ref().to_owned())))],
            idx: AtomicUsize::new(0),
            interleaved: 1,
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
            .collect();

        Ok(Self {
            readers,
            idx: AtomicUsize::new(0),
            interleaved: 1,
        })
    }

    /// Stream reads created from interleaved fastq records from an input file.
    pub fn from_file_interleaved(file: impl AsRef<str>, interleaved: usize) -> Result<Self> {
        let reader = Mutex::new(parse_fastx_file(file.as_ref()).map_err(|e| Error::FileIo {
            file: file.as_ref().to_owned(),
            source: Box::new(e),
        })?);

        Ok(Self {
            readers: vec![(reader, Arc::new(Origin::File(file.as_ref().to_owned())))],
            idx: AtomicUsize::new(0),
            interleaved,
        })
    }

    /// Stream reads created from fastq records from an arbitrary `Read`er.
    pub fn from_reader(reader: impl std::io::Read + Send + 'reader) -> Result<Self> {
        let reader =
            Mutex::new(parse_fastx_reader(reader).map_err(|e| Error::BytesIo(Box::new(e)))?);

        Ok(Self {
            readers: vec![(reader, Arc::new(Origin::Bytes))],
            idx: AtomicUsize::new(0),
            interleaved: 1,
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

        Ok(Self {
            readers,
            idx: AtomicUsize::new(0),
            interleaved: 1,
        })
    }

    /// Stream reads created from interleaved fastq records from an arbitrary `Read`er.
    pub fn from_interleaved_reader(
        reader: impl std::io::Read + Send + 'reader,
        interleaved: usize,
    ) -> Result<Self> {
        let reader =
            Mutex::new(parse_fastx_reader(reader).map_err(|e| Error::BytesIo(Box::new(e)))?);

        Ok(Self {
            readers: vec![(reader, Arc::new(Origin::Bytes))],
            idx: AtomicUsize::new(0),
            interleaved,
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
                        &record.seq(),
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
                        &record.seq(),
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
}
