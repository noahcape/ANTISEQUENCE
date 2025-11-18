use needletail::*;

use thread_local::*;

use std::cell::RefCell;
use std::collections::VecDeque;
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
            .unwrap_or(1024)
    })
}

pub struct InputFastqOp<'reader> {
    readers: Vec<(Mutex<Box<dyn FastxReader + 'reader>>, Arc<Origin>)>,
    buf: ThreadLocal<RefCell<VecDeque<Read>>>,
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
            buf: ThreadLocal::new(),
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
            buf: ThreadLocal::new(),
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
            buf: ThreadLocal::new(),
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
            buf: ThreadLocal::new(),
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
            buf: ThreadLocal::new(),
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
            buf: ThreadLocal::new(),
            idx: AtomicUsize::new(0),
            interleaved,
        })
    }
}

impl<'reader, T: Trace> GraphNode<T> for InputFastqOp<'reader> {
    fn run(&self, read: Option<Read>, trace: &T) -> Result<(Option<Read>, bool)> {
        let start = trace.start(&read);
        assert!(read.is_none(), "Expected no input reads for {}", Self::NAME);

        let cs = chunk_size();
        let buf = self
            .buf
            .get_or(|| RefCell::new(VecDeque::with_capacity(cs)));
        let mut b = buf.borrow_mut();

        if b.is_empty() {
            // Lock readers once per refill to amortize lock overhead
            let mut locked_readers = self
                .readers
                .iter()
                .map(|(r, o)| (r.lock(), o))
                .collect::<Vec<_>>();

            'outer: for _ in 0..cs {
                let idx = self.idx.fetch_add(self.interleaved, Ordering::Relaxed);
                let mut curr_read = Read::new();

                if self.interleaved > 1 {
                    // interleaved records all come from one file
                    let (locked_reader, origin) = &mut locked_readers[0];

                    for i in 0..self.interleaved {
                        let Some(record) = locked_reader.next() else {
                            if i == 0 {
                                break 'outer;
                            }
                            Err(Error::UnpairedRead(format!("\"{}\"", &**origin)))?
                        };
                        let record = record.map_err(|e| Error::ParseRecord {
                            origin: (***origin).clone(),
                            idx: idx + i,
                            source: Box::new(e),
                        })?;
                        let name_opt = Some(record.id());
                        let qual_opt = record.qual();
                        curr_read.add_fastq_parts(
                            (i + 1) as _,
                            name_opt,
                            &record.seq(),
                            qual_opt,
                            Arc::clone(origin),
                            idx + i,
                        );
                    }
                } else {
                    // gather records from multiple different files
                    for (i, (locked_reader, origin)) in locked_readers.iter_mut().enumerate() {
                        let Some(record) = locked_reader.next() else {
                            if i == 0 {
                                break 'outer;
                            }
                            Err(Error::UnpairedRead(format!("\"{}\"", &**origin)))?
                        };
                        let record = record.map_err(|e| Error::ParseRecord {
                            origin: (***origin).clone(),
                            idx,
                            source: Box::new(e),
                        })?;
                        let name_opt = Some(record.id());
                        let qual_opt = record.qual();
                        curr_read.add_fastq_parts(
                            (i + 1) as _,
                            name_opt,
                            &record.seq(),
                            qual_opt,
                            Arc::clone(origin),
                            idx,
                        );
                    }
                }

                b.push_back(curr_read);
            }
        }

        if b.is_empty() {
            return Ok((None, true));
        }

        let res = b.pop_front();
        trace.add(<Self as GraphNode<T>>::name(self), start, &res);
        Ok((res, false))
    }

    fn required_names(&self) -> &[LabelOrAttr] {
        &[]
    }

    fn name(&self) -> &'static str {
        Self::NAME
    }
}
