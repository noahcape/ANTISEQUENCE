use std::fs::File;
use std::io::{BufWriter, Write, IoSlice};
use std::sync::Arc;
use std::sync::OnceLock;
use std::borrow::Cow;
use parking_lot::Mutex;

use rustc_hash::FxHashMap;

use flate2::{write::GzEncoder, Compression};

use crate::graph::*;

struct BufState { buf: Vec<u8>, count: usize }

struct TlsOutputState {
    writers: FxHashMap<Vec<u8>, Arc<Mutex<dyn Write + Send>>>,
    bufs: FxHashMap<Vec<u8>, BufState>,
}

impl TlsOutputState {
    fn new() -> Self { Self { writers: FxHashMap::default(), bufs: FxHashMap::default() } }
}

impl Drop for TlsOutputState {
    fn drop(&mut self) {
        for (k, bs) in self.bufs.iter_mut() {
            if bs.buf.is_empty() { continue; }
            if let Some(w) = self.writers.get(k) {
                let mut w = w.lock();
                let _ = (&mut *w).write_all(&bs.buf);
                let _ = (&mut *w).flush();
                bs.buf.clear();
                bs.count = 0;
            }
        }
    }
}

thread_local! {
    static OUTPUT_TLS: std::cell::RefCell<TlsOutputState> = std::cell::RefCell::new(TlsOutputState::new());
}

pub struct OutputFastqFileOp {
    required_names: Vec<LabelOrAttr>,
    file_exprs: Vec<Expr>,
    file_consts: Vec<Option<Vec<u8>>>,
    file_writers: Mutex<FxHashMap<Vec<u8>, Arc<Mutex<dyn Write + Send>>>>,
}

impl OutputFastqFileOp {
    const NAME: &'static str = "OutputFastqFileOp";

    /// Output reads (read 1 only) to a file whose path is specified by an expression.
    pub fn from_file(file_expr: impl Into<Expr>) -> Self {
        let mut file_expr: Expr = file_expr.into();
        let _ = file_expr.optimize();

        let required_names = file_expr.required_names();

        let file_const = if required_names.is_empty() {
            let tmp_read = crate::read::Read::new();
            let const_bytes = file_expr
                .eval_bytes(&tmp_read, false)
                .unwrap_or_else(|e| panic!("{e}"))
                .into_owned();
            Some(const_bytes)
        } else {
            None
        };

        Self {
            required_names,
            file_exprs: vec![file_expr],
            file_consts: vec![file_const],
            file_writers: Mutex::new(FxHashMap::default()),
        }
    }

    /// Output reads to separate files whose paths are specified by expressions.
    pub fn from_files<E: Into<Expr>>(file_exprs: impl IntoIterator<Item = E>) -> Self {
        let mut file_exprs: Vec<Expr> = file_exprs.into_iter().map(|e| e.into()).collect::<Vec<_>>();

        for e in file_exprs.iter_mut() { let _ = e.optimize(); }
        let required_names = file_exprs
            .iter()
            .flat_map(|e| e.required_names().into_iter())
            .collect::<Vec<_>>();

        let file_consts = file_exprs
            .iter()
            .map(|e| {
                if e.required_names().is_empty() {
                    let tmp_read = crate::read::Read::new();
                    let const_bytes = e
                        .eval_bytes(&tmp_read, false)
                        .unwrap_or_else(|e| panic!("{e}"))
                        .into_owned();
                    Some(const_bytes)
                } else {
                    None
                }
            })
            .collect::<Vec<_>>();

        Self {
            required_names,
            file_exprs,
            file_consts,
            file_writers: Mutex::new(FxHashMap::default()),
        }
    }

    // get the corresponding file writer for each read first so writing to different files can be parallelized
    fn get_writer(&self, file_name: &[u8]) -> std::io::Result<Arc<Mutex<dyn Write + Send>>> {
        use std::collections::hash_map::Entry::*;
        let mut file_writers = self.file_writers.lock();

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

    fn get_cached_writer(&self, file_name: &[u8]) -> std::io::Result<Arc<Mutex<dyn Write + Send>>> {
        // if TLS is disabled, just create a new writer
        if writer_tls_disabled() { 
            return self.get_writer(file_name); 
        }
        // otherwise, try to get the writer from TLS
        if let Some(w) = OUTPUT_TLS.with(|m| m.borrow().writers.get(file_name).map(Arc::clone)) { return Ok(w); }
        // if not found, create a new writer and cache it
        let w = self.get_writer(file_name)?; 
        OUTPUT_TLS.with(|m| { 
            m.borrow_mut().writers.insert(file_name.to_vec(), Arc::clone(&w)); 
        });
        Ok(w)
    }
}

#[inline(always)]
fn stub_output() -> bool {
    static STUB: OnceLock<bool> = OnceLock::new();
    *STUB.get_or_init(|| {
        std::env::var("ANTISEQ_STUB_OUTPUT")
            .ok()
            .map(|v| v == "1" || v.eq_ignore_ascii_case("true"))
            .unwrap_or(false)
    })
}

#[inline(always)]
fn writer_tls_disabled() -> bool {
    static DIS: OnceLock<bool> = OnceLock::new();
    *DIS.get_or_init(|| {
        std::env::var("ANTISEQ_DISABLE_OUTPUT_TLS")
            .ok()
            .map(|v| v == "1" || v.eq_ignore_ascii_case("true"))
            .unwrap_or(false)
    })
}

#[inline(always)]
fn output_batch_size() -> Option<usize> {
    static BSZ: OnceLock<Option<usize>> = OnceLock::new();
    *BSZ.get_or_init(|| {
        std::env::var("ANTISEQ_OUTPUT_BATCH").ok()
            .and_then(|s| s.trim().parse::<usize>().ok())
            .and_then(|n| if n > 0 { Some(n) } else { None })
    })
}

impl<T: Trace> GraphNode<T> for OutputFastqFileOp {
    fn run_inner(&self, read: Read) -> Result<(Option<Read>, bool)> {
        if stub_output() { return Ok((Some(read), false)); }
        for (i, file_expr) in self.file_exprs.iter().enumerate() {
            let file_name: Cow<[u8]> = if let Some(c) = self.file_consts.get(i).and_then(|o| o.as_ref()) {
                Cow::Borrowed(&c[..])
            } else {
                file_expr
                    .eval_bytes(&read, false)
                    .map_err(|e| Error::NameError {
                        source: e,
                        read: read.clone(),
                        context: Self::NAME,
                    })?
            };

            let locked_writer = self.get_cached_writer(&file_name).map_err(|e| Error::FileIo { file: utf8(&file_name), source: Box::new(e) })?;

            let record = read.to_fastq((i + 1) as _).map_err(|e| Error::NameError { source: e, read: read.clone(), context: Self::NAME })?;

            if let Some(n) = output_batch_size() {
                let (name, seq, qual) = record;
                let mut tmp = Vec::with_capacity(1 + name.len() + 1 + seq.len() + 2 + 1 + qual.len() + 1);
                tmp.push(b'@'); tmp.extend_from_slice(name); tmp.push(b'\n');
                tmp.extend_from_slice(seq); tmp.push(b'\n'); tmp.push(b'+'); tmp.push(b'\n');
                tmp.extend_from_slice(qual); tmp.push(b'\n');

                let key = file_name.to_vec();
                let mut flush_needed = false;
                OUTPUT_TLS.with(|m| {
                    let mut s = m.borrow_mut();
                    let bs = s.bufs.entry(key.clone()).or_insert_with(|| BufState { buf: Vec::new(), count: 0 });
                    bs.buf.extend_from_slice(&tmp);
                    bs.count += 1;
                    if bs.count >= n { flush_needed = true; }
                });

                if flush_needed {
                    let mut to_write = Vec::new();
                    OUTPUT_TLS.with(|m| {
                        let mut s = m.borrow_mut();
                        if let Some(bs) = s.bufs.get_mut(&key) {
                            to_write = std::mem::take(&mut bs.buf);
                            bs.count = 0;
                        }
                    });
                    if !to_write.is_empty() {
                        let mut w = locked_writer.lock();
                        (&mut *w).write_all(&to_write).map_err(|e| Error::FileIo { file: utf8(&file_name), source: Box::new(e) })?;
                        (&mut *w).flush().map_err(|e| Error::FileIo { file: utf8(&file_name), source: Box::new(e) })?;
                    }
                }
            } else {
                let mut writer = locked_writer.lock();
                write_fastq_record(&mut *writer, record);
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
        if stub_output() { return Ok((Some(read), false)); }
        for (i, writer) in self.writers.iter().enumerate() {
            let record = read.to_fastq((i + 1) as _).map_err(|e| Error::NameError {
                source: e,
                read: read.clone(),
                context: Self::NAME,
            })?;

            let mut writer = writer.lock();
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

#[inline(always)]
pub fn write_fastq_record(
    writer: &mut (dyn Write + std::marker::Send),
    record: (&[u8], &[u8], &[u8]),
) {
    let (name, seq, qual) = record;
    let segs: [&[u8]; 7] = [b"@", name, b"\n", seq, b"\n+\n", qual, b"\n"];
    let total = 1 + name.len() + 1 + seq.len() + 3 + qual.len() + 1;

    let mut idx = 0usize; // current segment
    let mut off = 0usize; // offset into current segment
    let mut written = 0usize;

    while written < total {
        // Build IoSlices for remaining segments
        let mut buf: [IoSlice; 7] = [
            IoSlice::new(b""),
            IoSlice::new(b""),
            IoSlice::new(b""),
            IoSlice::new(b""),
            IoSlice::new(b""),
            IoSlice::new(b""),
            IoSlice::new(b""),
        ];
        let mut n = 0usize;
        let mut j = idx;
        while j < segs.len() {
            let s = segs[j];
            let slice = if j == idx { &s[off..] } else { s };
            if !slice.is_empty() {
                buf[n] = IoSlice::new(slice);
                n += 1;
            }
            j += 1;
        }

        match writer.write_vectored(&buf[..n]) {
            Ok(0) => {
                // Fallback: write some from current segment
                if idx >= segs.len() { break; }
                let first = &segs[idx][off..];
                if !first.is_empty() {
                    let nw = writer.write(first).unwrap();
                    if nw == 0 { continue; }
                    written += nw;
                    off += nw;
                    if off == segs[idx].len() { idx += 1; off = 0; }
                } else {
                    idx += 1; off = 0;
                }
            }
            Ok(nw) => {
                written += nw;
                let mut rem = nw;
                while rem > 0 {
                    let remain_in_cur = segs[idx].len() - off;
                    if rem < remain_in_cur {
                        off += rem;
                        rem = 0;
                    } else {
                        rem -= remain_in_cur;
                        idx += 1;
                        off = 0;
                    }
                }
            }
            Err(ref e) if e.kind() == std::io::ErrorKind::Interrupted => continue,
            Err(e) => panic!("{}", e),
        }
    }
}
