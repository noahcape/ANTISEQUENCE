use std::fs::File;
use std::io::{BufWriter, Write, IoSlice};
use std::sync::{Arc, Mutex};

use rustc_hash::FxHashMap;

use flate2::{write::GzEncoder, Compression};

use crate::graph::*;

pub struct OutputFastqFileOp {
    required_names: Vec<LabelOrAttr>,
    file_exprs: Vec<Expr>,
    file_writers: Mutex<FxHashMap<Vec<u8>, Arc<Mutex<dyn Write + Send>>>>,
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
                    let gz = GzEncoder::new(File::create(file_path)?, Compression::default());
                    Arc::new(Mutex::new(BufWriter::with_capacity(1 << 20, gz)))
                } else {
                    Arc::new(Mutex::new(BufWriter::with_capacity(1 << 20, File::create(file_path)?)))
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

            let locked_writer = self.get_writer(&file_name).map_err(|e| Error::FileIo {
                file: utf8(&file_name),
                source: Box::new(e),
            })?;

            let record = read.to_fastq((i + 1) as _).map_err(|e| Error::NameError {
                source: e,
                read: read.clone(),
                context: Self::NAME,
            })?;

            let mut writer = locked_writer.lock().unwrap();
            write_fastq_record(&mut *writer, record);
        }

        Ok((Some(read), false))
    }

    fn required_names(&self) -> &[LabelOrAttr] {
        &self.required_names
    }

    fn name(&self) -> &'static str {
        Self::NAME
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
}

pub fn write_fastq_record(
    writer: &mut (dyn Write + std::marker::Send),
    record: (&[u8], &[u8], &[u8]),
) {
    let (name, seq, qual) = record;

    let mut parts = [
        IoSlice::new(b"@" as &[u8]),
        IoSlice::new(name),
        IoSlice::new(b"\n" as &[u8]),
        IoSlice::new(seq),
        IoSlice::new(b"\n+\n" as &[u8]),
        IoSlice::new(qual),
        IoSlice::new(b"\n" as &[u8]),
    ];

    let total: usize = parts.iter().map(|p| p.len()).sum();
    let n = writer.write_vectored(&parts).unwrap();
    if n == total {
        return;
    }
    if n == 0 {
        // Fallback to ensure forward progress
        for p in &parts {
            writer.write_all(p.as_ref()).unwrap();
        }
        return;
    }

    // Write remaining bytes sequentially
    let mut remaining = n;
    for (idx, p) in parts.iter().enumerate() {
        let s = p.as_ref();
        if remaining < s.len() {
            writer.write_all(&s[remaining..]).unwrap();
            for k in idx + 1..parts.len() {
                writer.write_all(parts[k].as_ref()).unwrap();
            }
            return;
        } else {
            remaining -= s.len();
        }
    }
}
