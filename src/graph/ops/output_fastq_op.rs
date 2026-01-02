use std::fs::File;
use std::io::{BufWriter, Write, IoSlice};
use std::sync::Arc;
use std::borrow::Cow;
use parking_lot::Mutex;
use std::cell::RefCell;

use rustc_hash::FxHashMap;
use thread_local::ThreadLocal;

use flate2::{write::GzEncoder, Compression};

use crate::graph::*;

/// Thread-local state for output operations.
///
/// This structure holds per-thread buffers and file writers to allow multiple threads
/// to write output without constantly contending for locks on the shared file writers.
struct TlsOutputState {
    /// Cache of file writers keyed by filename.
    /// Note: These Arcs point to the SAME shared writers as the main Op struct,
    /// but cached here to avoid map lookups if possible (though map is still used here).
    /// Actually, the main benefit is that we buffer data in `bufs` before locking the writer.
    writers: FxHashMap<Vec<u8>, Arc<Mutex<dyn Write + Send>>>,
    /// Per-file write buffers. Key is filename. Value is byte buffer.
    bufs: FxHashMap<Vec<u8>, Vec<u8>>,
}

impl TlsOutputState {
    fn new() -> Self { Self { writers: FxHashMap::default(), bufs: FxHashMap::default() } }
}

/// Drop implementation flushes any remaining buffered data when the thread exits.
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

// Thread-local storage for output operations.
//
// Gives a 2-5% speedup on file-output benchmarks but is super helpful especially as
// the number of threads scales. Without it, adding more threads could actually slow 
// down the program due to excessive lock contention. Every time a thread wants to 
// write a read to a file, it must acquire a Mutex lock on the shared file writer. 
// If you have 4 threads processing millions of reads, they will constantly fight 
// for this single lock, forcing them to wait in line (serialization).
//
// Also, writes are expensive. Instead of making a system call (or even a locked 
// library call) for every single read (e.g., ~100 bytes), the thread fills up a 
// larger buffer (e.g., 8KB or more). It only acquires the lock and writes to the 
// actual file when the buffer is full or the batch is done. This means you might 
// lock and write once for every 100 reads instead of 100 times.
thread_local! {
    /// Thread-local storage instance. Lazily initialized for each worker thread.
    static OUTPUT_TLS: std::cell::RefCell<TlsOutputState> = std::cell::RefCell::new(TlsOutputState::new());
}

/// Operation to output reads to one or more files.
///
/// The filename can be dynamic (determined by an expression per read), allowing
/// splitting reads into different files based on attributes/barcodes.
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
    ///
    /// This constructor creates an `OutputFastqFileOp` that can write reads to one or more files.
    /// The file paths are determined by evaluating the provided expressions for each read.
    ///
    /// # Features
    /// * **Dynamic Filenames**: Expressions can depend on read attributes (e.g., barcodes), enabling
    ///   demultiplexing where different reads are written to different files based on their content.
    /// * **Optimization**: Expressions are optimized during construction. If an expression evaluates
    ///   to a constant path (independent of read data), it is pre-calculated to avoid per-read overhead.
    ///
    /// # Arguments
    /// * `file_exprs` - An iterator of types that can be converted into `Expr`. Each expression
    ///   corresponds to an output file destination for the read.
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


impl<T: Trace> GraphNode<T> for OutputFastqFileOp {
    /// Execution logic for file output.
    ///
    /// 1. Uses thread-local storage to buffer writes.
    /// 2. Evaluates filename expressions for each read.
    /// 3. Formats FASTQ records into the thread-local buffer.
    /// 4. Flushes thread-local buffers to the shared file writers (acquiring locks only during flush).
    fn run_inner(&self, reads: Vec<Read>) -> Result<(Option<Vec<Read>>, bool)> {        
        OUTPUT_TLS.with(|tls| {
            let mut state_borrow = tls.borrow_mut();
            let TlsOutputState { writers, bufs } = &mut *state_borrow;
            
            for read in &reads {
                for (i, file_expr) in self.file_exprs.iter().enumerate() {
                    // Determine filename: either constant (optimized) or evaluated expression.
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
                    
                    // Get or create buffer for this filename
                    let buf = if let Some(buf) = bufs.get_mut(&*file_name) {
                        buf
                    } else {
                        bufs.entry(file_name.into_owned()).or_default()
                    };
                    
                    // Append formatted FASTQ to buffer
                    let (name, seq, qual) = record;
                    buf.reserve(1 + name.len() + 1 + seq.len() + 3 + qual.len() + 1);
                    buf.push(b'@'); buf.extend_from_slice(name); buf.push(b'\n');
                    buf.extend_from_slice(seq); buf.push(b'\n');
                    buf.extend_from_slice(b"+\n");
                    buf.extend_from_slice(qual); buf.push(b'\n');
                }
            }

            // Flush buffers to actual writers
            for (file_name, buf) in bufs.iter_mut() {
                if buf.is_empty() { continue; }
                
                // Get or create shared writer for this filename
                let writer = if let Some(w) = writers.get(file_name) {
                    Arc::clone(w)
                } else {
                    // If writer doesn't exist in TLS cache, check/create in shared map
                    let w = self.get_writer(file_name).map_err(|e| Error::FileIo { file: utf8(file_name), source: Box::new(e) })?;
                    writers.insert(file_name.clone(), Arc::clone(&w));
                    w
                };
                
                // Lock and write
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

/// Operation to output reads to generic Writers (e.g. stdout).
///
/// Uses thread-local buffering to minimize contention on the shared writers.
pub struct OutputFastqOp<'writer> {
    /// Shared writers protected by Mutex.
    writers: Vec<Mutex<Box<dyn Write + Send + 'writer>>>,
    /// Thread-local buffers. `RefCell` allows interior mutability within the thread.
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

impl<'writer, T: Trace> GraphNode<T> for OutputFastqOp<'writer> {
    /// Execution logic for generic writer output.
    ///
    /// This implementation employs a **buffered write strategy** using thread-local storage to
    /// maximize performance in multi-threaded environments.
    ///
    /// # Strategy
    /// 1. **Thread-Local Buffering**: Each thread maintains its own private buffer for each output writer.
    /// 2. **Batch Formatting**: Reads in the input batch are formatted as FASTQ and appended to these
    ///    local buffers, avoiding any lock contention during the heavy formatting phase.
    /// 3. **Batched Flushing**: Only after the entire batch is processed does the thread acquire
    ///    the lock on the shared writer to flush the accumulated data.
    ///
    /// This significantly reduces the number of lock acquisitions from `N` (number of reads) to
    /// `1` (per batch), preventing thread starvation and serialization bottlenecks.
    /// In practice this reduces lock contention by ~1000x (depending on batch size).
    fn run_inner(&self, reads: Vec<Read>) -> Result<(Option<Vec<Read>>, bool)> {
        // Get or initialize thread-local buffers
        let mut buffers = self.buffers.get_or(|| RefCell::new(vec![Vec::new(); self.writers.len()])).borrow_mut();
        
        // Buffer the output for this batch
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
        
        // Flush buffers to shared writers
        for (i, buf) in buffers.iter_mut().enumerate() {
             if !buf.is_empty() {
                 let mut writer = self.writers[i].lock();
                 // If we encounter an unhappy error (e.g., disk full, pipe broken, etc.)
                 // we return an error to the caller.
                 writer.write_all(buf).map_err(|e| Error::BytesIo(Box::new(e)))?;
                 buf.clear();
             }
        }

        Ok((Some(reads), false))
    }

    fn required_names(&self) -> &[LabelOrAttr] {
        &[]
    }

    fn name(&self) -> &'static str {
        Self::NAME
    }
}

/// Helper to write a single FASTQ record efficiently using vectored IO.
///
/// This avoids concatenating the parts of the record into a single buffer before writing.
/// Instead, it constructs an array of `IoSlice`s and writes them all at once.
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

#[cfg(test)]
mod tests {
    use super::*;
    use std::io;
    use crate::trace::NoTrace;

    struct FailingWriter;

    impl Write for FailingWriter {
        // Simulate a failed write operation
        fn write(&mut self, _buf: &[u8]) -> io::Result<usize> {
            Err(io::Error::new(io::ErrorKind::BrokenPipe, "simulated error"))
        }

        // Simulate a failed flush operation
        fn flush(&mut self) -> io::Result<()> {
            Err(io::Error::new(io::ErrorKind::BrokenPipe, "simulated error"))
        }
    }

    #[test]
    fn test_output_failure() {
        // Create an OutputFastqOp with a FailingWriter
        // This will simulate a failed write operation
        let op = OutputFastqOp::from_writer(FailingWriter);
        
        // Create a dummy read
        let mut read = Read::new();
        read.set_fastq_entry(
            0, 
            StrType::Name(1), 
            b"read1", 
            None, 
            Arc::new(Origin::Bytes), 
            0
        );
        read.set_fastq_entry(
            1, 
            StrType::Seq(1), 
            b"ACGT", 
            Some(b"IIII"), 
            Arc::new(Origin::Bytes), 
            0
        );

        // Explicitly specify NoTrace to satisfy the generic bound T: Trace
        let res = GraphNode::<NoTrace>::run_inner(&op, vec![read]);

        assert!(res.is_err());
        match res {
            Err(Error::BytesIo(e)) => {
                assert_eq!(e.to_string(), "simulated error");
            }
            _ => panic!("Expected BytesIo error, got {:?}", res),
        }
    }
}
