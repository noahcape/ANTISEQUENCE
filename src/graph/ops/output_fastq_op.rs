use std::fs::File;
use std::io::{BufWriter, Write, IoSlice};
use std::sync::Arc;
use std::sync::OnceLock;
use std::borrow::Cow;
use parking_lot::Mutex;
use std::cell::RefCell;

use rustc_hash::FxHashMap;
use thread_local::ThreadLocal;

use flate2::{write::GzEncoder, Compression};

use crate::graph::*;

struct TlsOutputState {
    writers: FxHashMap<Vec<u8>, Arc<Mutex<dyn Write + Send>>>,
    bufs: FxHashMap<Vec<u8>, Vec<u8>>,
}

impl TlsOutputState {
    fn new() -> Self { Self { writers: FxHashMap::default(), bufs: FxHashMap::default() } }
}

impl Drop for TlsOutputState {
    fn drop(&mut self) {
        for (k, buf) in self.bufs.iter_mut() {
            if buf.is_empty() { continue; }
            if let Some(w) = self.writers.get(k) {
                let mut w = w.lock();
                let _ = (&mut *w).write_all(buf);
                buf.clear();
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
}

impl Drop for OutputFastqFileOp {
    fn drop(&mut self) {
        let writers = self.file_writers.lock();
        for writer in writers.values() {
            let mut w = writer.lock();
            let _ = w.flush();
        }
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

impl<T: Trace> GraphNode<T> for OutputFastqFileOp {
    fn run_inner(&self, reads: Vec<Read>) -> Result<(Option<Vec<Read>>, bool)> {
        if stub_output() { return Ok((Some(reads), false)); }
        
        OUTPUT_TLS.with(|tls| {
            let mut state_borrow = tls.borrow_mut();
            let TlsOutputState { writers, bufs } = &mut *state_borrow;
            
            for read in &reads {
                for (i, file_expr) in self.file_exprs.iter().enumerate() {
                    let file_name: Cow<[u8]> = if let Some(c) = self.file_consts.get(i).and_then(|o| o.as_ref()) {
                        Cow::Borrowed(&c[..])
                    } else {
                        file_expr
                            .eval_bytes(read, false)
                            .map_err(|e| Error::NameError {
                                source: e,
                                read: read.clone(),
                                context: Self::NAME,
                            })?
                    };

                    let record = read.to_fastq((i + 1) as _).map_err(|e| Error::NameError { source: e, read: read.clone(), context: Self::NAME })?;
                    
                    let buf = if let Some(buf) = bufs.get_mut(&*file_name) {
                        buf
                    } else {
                        bufs.entry(file_name.into_owned()).or_default()
                    };
                    
                    let (name, seq, qual) = record;
                    buf.reserve(1 + name.len() + 1 + seq.len() + 3 + qual.len() + 1);
                    buf.push(b'@'); buf.extend_from_slice(name); buf.push(b'\n');
                    buf.extend_from_slice(seq); buf.push(b'\n');
                    buf.extend_from_slice(b"+\n");
                    buf.extend_from_slice(qual); buf.push(b'\n');
                }
            }

            // Flush buffers
            for (file_name, buf) in bufs.iter_mut() {
                if buf.is_empty() { continue; }
                
                let writer = if let Some(w) = writers.get(file_name) {
                    Arc::clone(w)
                } else {
                    let w = self.get_writer(file_name).map_err(|e| Error::FileIo { file: utf8(file_name), source: Box::new(e) })?;
                    writers.insert(file_name.clone(), Arc::clone(&w));
                    w
                };
                
                let mut w = writer.lock();
                w.write_all(buf).map_err(|e| Error::FileIo { file: utf8(file_name), source: Box::new(e) })?;
                
                buf.clear();
            }
            
            Ok::<(), Error>(())
        }).map_err(|e| e)?; // Extract result from with() which returns whatever closure returns.
        // Wait, `with` returns R. My closure returns Result<(), Error>.
        // So map_err is correct if I propagate it.

        Ok((Some(reads), false))
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
    buffers: ThreadLocal<RefCell<Vec<Vec<u8>>>>,
}

impl<'writer> OutputFastqOp<'writer> {
    const NAME: &'static str = "OutputFastqOp";

    /// Output reads (read 1 only) to a `Write`r.
    pub fn from_writer(writer: impl Write + Send + 'writer) -> Self {
        Self {
            writers: vec![Mutex::new(Box::new(writer))],
            buffers: ThreadLocal::new(),
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
            buffers: ThreadLocal::new(),
        }
    }
}

impl<'writer> Drop for OutputFastqOp<'writer> {
    fn drop(&mut self) {
        // Explicitly flush all writers to ensure data is written before file close.
        // BufWriter::drop can silently ignore errors, so we flush explicitly here.
        for writer in &self.writers {
            let mut w = writer.lock();
            let _ = w.flush();
        }
    }
}

impl<'writer, T: Trace> GraphNode<T> for OutputFastqOp<'writer> {
    fn run_inner(&self, reads: Vec<Read>) -> Result<(Option<Vec<Read>>, bool)> {
        if stub_output() { return Ok((Some(reads), false)); }
        
        let mut buffers = self.buffers.get_or(|| RefCell::new(vec![Vec::new(); self.writers.len()])).borrow_mut();
        
        for read in &reads {
            for (i, buf) in buffers.iter_mut().enumerate() {
                let record = read.to_fastq((i + 1) as _).map_err(|e| Error::NameError {
                    source: e,
                    read: read.clone(),
                    context: Self::NAME,
                })?;
                
                let (name, seq, qual) = record;
                buf.reserve(1 + name.len() + 1 + seq.len() + 3 + qual.len() + 1);
                buf.push(b'@'); buf.extend_from_slice(name); buf.push(b'\n');
                buf.extend_from_slice(seq); buf.push(b'\n');
                buf.extend_from_slice(b"+\n");
                buf.extend_from_slice(qual); buf.push(b'\n');
            }
        }
        
        // Lock ALL writers at once to ensure R1 and R2 are written atomically
        // This prevents interleaving issues with multi-threaded output
        let mut locked_writers: Vec<_> = self.writers.iter().map(|w| w.lock()).collect();
        
        for (i, buf) in buffers.iter_mut().enumerate() {
             if !buf.is_empty() {
                 locked_writers[i].write_all(buf).unwrap(); // TODO: proper error handling
                 buf.clear();
             }
        }
        // All locks released together when locked_writers is dropped

        Ok((Some(reads), false))
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
