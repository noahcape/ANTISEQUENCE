use std::{
    fs::File,
    io::{BufWriter, Write},
    path::Path,
    sync::Mutex,
    time::{Duration, Instant},
};

use serde::Serialize;
use serde_json;

use crate::read::*;

pub static DEFAULT_TRACE_PATH: &'static str = "ANTISEQUENCE.trace.json";

pub trait Trace: Send + Sync {
    type S;

    fn new(file_path: impl AsRef<Path>) -> Self;
    fn start(&self, reads: &Option<Vec<Read>>) -> Self::S;
    fn add(&self, name: &str, start: Self::S, reads: &Option<Vec<Read>>);
    fn finish(self);
}

pub struct TraceReads {
    start: Instant,
    writer: Mutex<(bool, BufWriter<File>)>,
}

impl Trace for TraceReads {
    type S = (Duration, Instant, Option<usize>);

    fn new(file_path: impl AsRef<Path>) -> Self {
        let path = file_path.as_ref();
        if let Some(parent) = path.parent() {
            std::fs::create_dir_all(parent).unwrap();
        }

        let mut writer = BufWriter::with_capacity(1 << 20, File::create(path).unwrap());
        writer.write_all(BEFORE).unwrap();
        let writer = Mutex::new((true, writer));

        Self {
            start: Instant::now(),
            writer,
        }
    }

    fn start(&self, reads: &Option<Vec<Read>>) -> Self::S {
        (
            self.start.elapsed(),
            Instant::now(),
            reads.as_ref().and_then(|v| v.first()).map(|r| r.first_idx()),
        )
    }

    fn add(&self, name: &str, start: Self::S, reads: &Option<Vec<Read>>) {
        let start_us = (start.0.as_nanos() as f64) / 1000.0f64;
        let dur_us = (start.1.elapsed().as_nanos() as f64) / 1000.0f64;
        let first_idx = reads.as_ref().and_then(|v| v.first()).map(|r| r.first_idx()).or(start.2).unwrap();
        let event = TraceEvent::new(name, start_us, dur_us, first_idx, reads);
        let mut writer = self.writer.lock().unwrap();

        if !writer.0 {
            writer.1.write_all(b",\n").unwrap();
        }

        writer.0 = false;
        serde_json::to_writer(&mut writer.1, &event).unwrap();
    }

    fn finish(self) {
        let mut writer = self.writer.into_inner().unwrap();
        writer.1.write_all(AFTER).unwrap();
    }
}

pub struct NoTrace;

impl Trace for NoTrace {
    type S = ();

    fn new(_file_path: impl AsRef<Path>) -> Self {
        Self
    }

    fn start(&self, _reads: &Option<Vec<Read>>) -> Self::S {
        ()
    }

    fn add(&self, _name: &str, _start: Self::S, _reads: &Option<Vec<Read>>) {}
    fn finish(self) {}
}

static BEFORE: &'static [u8] = br#"{
  "traceEvents": [
"#;

static AFTER: &'static [u8] = br#"
  ]
}
"#;

#[derive(Serialize)]
struct TraceEvent<'a> {
    name: &'a str,
    ph: char,
    ts: f64,
    dur: f64,
    pid: usize,
    tid: usize,
    args: Args,
}

#[derive(Serialize)]
struct Args {
    read: Option<SerializableRead>,
}

impl<'a> TraceEvent<'a> {
    pub fn new(name: &'a str, start: f64, dur: f64, tid: usize, reads: &'a Option<Vec<Read>>) -> Self {
        Self {
            name,
            ph: 'X',
            ts: start,
            dur,
            pid: 0,
            tid,
            args: Args {
                read: reads.as_ref().and_then(|v| v.first()).map(|r| SerializableRead::from(r)),
            },
        }
    }
}
