use std::fs::File;
use std::io::{BufWriter, Write};
use std::sync::Mutex;

use serde_json;

use flate2::{write::GzEncoder, Compression};

use crate::graph::*;

pub struct OutputJsonNode<'writer> {
    writer: Mutex<Box<dyn Write + Send + 'writer>>,
}

impl<'writer> OutputJsonNode<'writer> {
    const NAME: &'static str = "OutputJsonNode";

    /// Output reads to a file in JSONL format.
    pub fn from_file(file: impl AsRef<str>) -> std::io::Result<Self> {
        let file_path = file.as_ref();

        if let Some(parent) = std::path::Path::new(file_path).parent() {
            std::fs::create_dir_all(parent)?;
        }

        let writer: Mutex<Box<dyn Write + Send>> = if file_path.ends_with(".gz") {
            Mutex::new(Box::new(BufWriter::new(GzEncoder::new(
                File::create(file_path)?,
                Compression::default(),
            ))))
        } else {
            Mutex::new(Box::new(BufWriter::new(File::create(file_path)?)))
        };

        Ok(Self { writer })
    }

    /// Output reads to a `Write`r in JSONL format.
    pub fn from_writer(writer: impl Write + Send + 'writer) -> Self {
        Self {
            writer: Mutex::new(Box::new(writer)),
        }
    }
}

impl<'writer> GraphNode for OutputJsonNode<'writer> {
    fn run(&self, read: Option<Read>) -> Result<(Option<Read>, bool)> {
        let Some(read) = read else {
            panic!("Expected some read!")
        };

        let mut writer = self.writer.lock().unwrap();
        serde_json::to_writer(&mut *writer, &SerializableRead::from(&read))
            .map_err(|e| Error::BytesIo(Box::new(e)))?;
        writeln!(&mut *writer).map_err(|e| Error::BytesIo(Box::new(e)))?;

        Ok((Some(read), false))
    }

    fn required_names(&self) -> &[LabelOrAttr] {
        &[]
    }

    fn name(&self) -> &'static str {
        Self::NAME
    }
}
