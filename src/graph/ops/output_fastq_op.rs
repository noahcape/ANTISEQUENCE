use std::cmp::Ordering;
use std::collections::hash_map::Entry;
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::sync::{Arc, Mutex};

use rustc_hash::FxHashMap;

use flate2::{write::GzEncoder, Compression};

use crate::graph::*;

const MEGABYTE: usize = 1000000;

pub struct OutputFastqFileOp {
    required_names: Vec<LabelOrAttr>,
    file_exprs: Vec<Expr>,
    file_writers: Mutex<FxHashMap<Vec<u8>, Arc<Mutex<dyn Write + Send>>>>,
    buffer: Mutex<HashMap<Vec<u8>, Vec<(Vec<u8>, Vec<u8>, Vec<u8>)>>>,
    buffer_size: Mutex<usize>,
}

impl OutputFastqFileOp {
    const NAME: &'static str = "OutputFastqFileOp";

    /// Output reads (read 1 only) to a file whose path is specified by an expression.
    pub fn from_file(file_expr: impl Into<Expr>) -> Self {
        let file_expr = file_expr.into();

        Self {
            required_names: file_expr.required_names(),
            file_exprs: vec![file_expr],
            file_writers: Mutex::new(FxHashMap::default()),
            buffer: Mutex::new(HashMap::new()),
            buffer_size: Mutex::new(0),
        }
    }

    /// Output reads to separate files whose paths are specified by expressions.
    pub fn from_files<E: Into<Expr>>(file_exprs: impl IntoIterator<Item = E>) -> Self {
        let file_exprs = file_exprs.into_iter().map(|e| e.into()).collect::<Vec<_>>();
        let required_names = file_exprs
            .iter()
            .flat_map(|e| e.required_names().into_iter())
            .collect::<Vec<_>>();

        Self {
            required_names,
            file_exprs,
            file_writers: Mutex::new(FxHashMap::default()),
            buffer: Mutex::new(HashMap::new()),
            buffer_size: Mutex::new(0),
        }
    }

    // get the corresponding file writer for each read first so writing to different files can be parallelized
    fn get_writer(&self, file_name: &[u8]) -> std::io::Result<Arc<Mutex<dyn Write + Send>>> {
        use std::collections::hash_map::Entry::*;
        let mut file_writers = self.file_writers.lock().unwrap();

        match file_writers.entry(file_name.to_owned()) {
            Occupied(e) => Ok(Arc::clone(e.get())),
            Vacant(e) => {
                // need to create the output file
                let file_path = std::str::from_utf8(file_name).unwrap();

                if let Some(parent) = std::path::Path::new(file_path).parent() {
                    std::fs::create_dir_all(parent)?;
                }

                let writer: Arc<Mutex<dyn Write + Send>> = if file_path.ends_with(".gz") {
                    Arc::new(Mutex::new(BufWriter::new(GzEncoder::new(
                        File::create(file_path)?,
                        Compression::default(),
                    ))))
                } else {
                    Arc::new(Mutex::new(BufWriter::new(File::create(file_path)?)))
                };

                Ok(Arc::clone(e.insert(writer)))
            }
        }
    }
}

impl<T: Trace> GraphNode<T> for OutputFastqFileOp {
    fn run_inner(&self, read: Read) -> Result<(Option<Read>, bool)> {
        for (i, file_expr) in self.file_exprs.iter().enumerate() {
            let file_name = file_expr
                .eval_bytes(&read, false)
                .map_err(|e| Error::NameError {
                    source: e,
                    read: read.clone(),
                    context: Self::NAME,
                })?;

            let record = read.to_fastq((i + 1) as _).map_err(|e| Error::NameError {
                source: e,
                read: read.clone(),
                context: Self::NAME,
            })?;

            let tupled_record = (record.0.to_vec(), record.1.to_vec(), record.2.to_vec());

            let mut locked_buffer = self.buffer.lock().unwrap();
            let mut locked_buffer_size = self.buffer_size.lock().unwrap();

            match locked_buffer.entry(file_name.to_vec()) {
                Entry::Occupied(mut e) => {
                    e.get_mut().push(tupled_record);
                }
                Entry::Vacant(e) => {
                    e.insert(vec![tupled_record]);
                }
            };

            *locked_buffer_size += record_size(record);

            if matches!(
                locked_buffer_size.cmp(&MEGABYTE),
                Ordering::Greater | Ordering::Equal
            ) {
                for fname in locked_buffer.keys() {
                    let locked_writer = self.get_writer(&fname).map_err(|e| Error::FileIo {
                        file: utf8(&fname),
                        source: Box::new(e),
                    })?;

                    let mut writer = locked_writer.lock().unwrap();
                    for (s, r, c) in locked_buffer.get(fname).unwrap() {
                        write_fastq_record(&mut *writer, (&s, &r, &c));
                    }
                }

                locked_buffer.clear();
                *locked_buffer_size = 0;
            }
        }

        Ok((Some(read), false))
    }

    fn required_names(&self) -> &[LabelOrAttr] {
        &self.required_names
    }

    fn name(&self) -> &'static str {
        Self::NAME
    }

    fn finish(&self) -> Result<bool> {
        let mut locked_buffer = self.buffer.lock().unwrap();
        for fname in locked_buffer.keys() {
            let locked_writer = self.get_writer(&fname).map_err(|e| Error::FileIo {
                file: utf8(&fname),
                source: Box::new(e),
            })?;

            let mut writer = locked_writer.lock().unwrap();
            for (s, r, c) in locked_buffer.get(fname).unwrap() {
                write_fastq_record(&mut *writer, (&s, &r, &c));
            }
        }

        locked_buffer.clear();
        *self.buffer_size.lock().unwrap() = 0;

        Ok(true)
    }
}

pub struct OutputFastqOp<'writer> {
    writers: Vec<Mutex<Box<dyn Write + Send + 'writer>>>,
}

impl<'writer> OutputFastqOp<'writer> {
    const NAME: &'static str = "OutputFastqOp";

    /// Output reads (read 1 only) to a `Write`r.
    pub fn from_writer(writer: impl Write + Send + 'writer) -> Self {
        Self {
            writers: vec![Mutex::new(Box::new(writer))],
        }
    }

    /// Output reads to separate `Write`rs.
    pub fn from_writers<W: Write + Send + 'writer>(writers: impl IntoIterator<Item = W>) -> Self {
        Self {
            writers: writers
                .into_iter()
                .map(|w| {
                    let w: Box<dyn Write + Send + 'writer> = Box::new(w);
                    Mutex::new(w)
                })
                .collect(),
        }
    }
}

impl<'writer, T: Trace> GraphNode<T> for OutputFastqOp<'writer> {
    fn run_inner(&self, read: Read) -> Result<(Option<Read>, bool)> {
        for (i, writer) in self.writers.iter().enumerate() {
            let record = read.to_fastq((i + 1) as _).map_err(|e| Error::NameError {
                source: e,
                read: read.clone(),
                context: Self::NAME,
            })?;

            let mut writer = writer.lock().unwrap();
            write_fastq_record(&mut *writer, record);
        }

        Ok((Some(read), false))
    }

    fn required_names(&self) -> &[LabelOrAttr] {
        &[]
    }

    fn name(&self) -> &'static str {
        Self::NAME
    }

    fn finish(&self) -> Result<bool> {
        Ok(true)
    }
}

pub fn write_fastq_record(
    writer: &mut (dyn Write + std::marker::Send),
    record: (&[u8], &[u8], &[u8]),
) {
    writer.write_all(b"@").unwrap();
    writer.write_all(&record.0).unwrap();
    writer.write_all(b"\n").unwrap();
    writer.write_all(&record.1).unwrap();
    writer.write_all(b"\n+\n").unwrap();
    writer.write_all(&record.2).unwrap();
    writer.write_all(b"\n").unwrap();
}

pub fn record_size(record: (&[u8], &[u8], &[u8])) -> usize {
    let (source, read, context) = record;
    core::mem::size_of_val(source) + core::mem::size_of_val(read) + core::mem::size_of_val(context)
}
