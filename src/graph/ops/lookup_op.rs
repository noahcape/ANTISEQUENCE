//! Lookup operation for mapping barcode values to sample identifiers.

use std::path::Path;
use std::fs::File;
use std::io::{BufRead, BufReader};

use rustc_hash::FxHashMap;

use crate::graph::*;
use crate::inline_string::InlineString;
use crate::read::Data;

/// Looks up a sequence value in a mapping table and sets an output attribute.
///
/// This is primarily used for demultiplexing, where a barcode sequence is looked up
/// in a sample mapping table to determine which sample the read belongs to.
pub struct LookupOp {
    required_names: Vec<LabelOrAttr>,
    /// The label to read the input sequence from (e.g., seq2.bc1)
    input_label: Label,
    /// The attribute name to set with the lookup result (e.g., "sample")
    output_attr: InlineString,
    /// Lookup table: barcode bytes -> sample name bytes
    lookup_table: FxHashMap<Vec<u8>, Vec<u8>>,
    /// Default value when barcode not found in table
    default_value: Vec<u8>,
}

impl LookupOp {
    const NAME: &'static str = "LookupOp";

    /// Create a new LookupOp.
    ///
    /// # Arguments
    /// * `input_label` - Label to read input sequence from (e.g., "seq2.bc1")
    /// * `output_attr` - Attribute name to store the lookup result (e.g., "sample")
    /// * `lookup_table` - HashMap mapping input bytes to output bytes
    /// * `default_value` - Value to use when input is not found in table
    pub fn new(
        input_label: impl Into<Label>,
        output_attr: &str,
        lookup_table: FxHashMap<Vec<u8>, Vec<u8>>,
        default_value: impl Into<Vec<u8>>,
    ) -> Self {
        let input_label = input_label.into();
        let required_names = vec![LabelOrAttr::Label(input_label.clone())];

        Self {
            required_names,
            input_label,
            output_attr: InlineString::new(output_attr.as_bytes()),
            lookup_table,
            default_value: default_value.into(),
        }
    }

    /// Load a lookup table from a TSV file.
    ///
    /// File format: two tab-separated columns per line.
    /// First column is the key (e.g., barcode), second is the value (e.g., sample name).
    /// Lines starting with '#' are treated as comments.
    pub fn load_tsv(path: impl AsRef<Path>) -> Result<FxHashMap<Vec<u8>, Vec<u8>>> {
        let file = File::open(path.as_ref()).map_err(|e| Error::FileIo {
            file: path.as_ref().to_string_lossy().to_string(),
            source: Box::new(e),
        })?;
        let reader = BufReader::new(file);
        let mut table = FxHashMap::default();

        for line in reader.lines() {
            let line = line.map_err(|e| Error::FileIo {
                file: path.as_ref().to_string_lossy().to_string(),
                source: Box::new(e),
            })?;
            let line = line.trim();

            // Skip empty lines and comments
            if line.is_empty() || line.starts_with('#') {
                continue;
            }

            let parts: Vec<&str> = line.split('\t').collect();
            if parts.len() >= 2 {
                table.insert(parts[0].as_bytes().to_vec(), parts[1].as_bytes().to_vec());
            }
        }

        Ok(table)
    }

    /// Create a LookupOp from a TSV file.
    pub fn from_tsv(
        input_label: impl Into<Label>,
        output_attr: &str,
        tsv_path: impl AsRef<Path>,
        default_value: impl Into<Vec<u8>>,
    ) -> Result<Self> {
        let lookup_table = Self::load_tsv(tsv_path)?;
        Ok(Self::new(input_label, output_attr, lookup_table, default_value))
    }
}

impl<T: Trace> GraphNode<T> for LookupOp {
    fn run_inner(&self, mut reads: Vec<Read>) -> Result<(Option<Vec<Read>>, bool)> {
        for read in &mut reads {
            // Get the input sequence bytes and clone to avoid borrow issues
            let input_bytes = match read.substring(self.input_label.str_type, self.input_label.label) {
                Ok(bytes) => bytes.to_vec(),
                Err(e) => {
                    return Err(Error::NameError {
                        source: e,
                        read: read.clone(),
                        context: Self::NAME,
                    });
                }
            };

            // Look up the value in the table
            let output_value = self
                .lookup_table
                .get(&input_bytes)
                .cloned()
                .unwrap_or_else(|| self.default_value.clone());

            // Set the output attribute
            match read.data_mut(self.input_label.str_type, self.input_label.label, self.output_attr) {
                Ok(data) => *data = Data::Bytes(output_value),
                Err(e) => {
                    return Err(Error::NameError {
                        source: e,
                        read: read.clone(),
                        context: Self::NAME,
                    });
                }
            }
        }

        Ok((Some(reads), false))
    }

    fn required_names(&self) -> &[LabelOrAttr] {
        &self.required_names
    }

    fn name(&self) -> &'static str {
        Self::NAME
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::read::{Origin, Read, StrType};
    use crate::trace::NoTrace;
    use std::sync::Arc;

    fn make_test_read(seq: &[u8]) -> Read {
        let mut read = Read::new();
        let origin = Arc::new(Origin::File("test.fastq".into()));
        read.add_fastq(2, b"test_read", seq, &vec![b'I'; seq.len()], origin, 0);
        
        // Add a mapping for bc1 label at position 0-8
        read.str_mappings_mut(StrType::Seq(2))
            .unwrap()
            .add_mapping(Some(InlineString::new(b"bc1")), 0, 8.min(seq.len()));
        
        read
    }

    #[test]
    fn test_lookup_found() {
        let mut table = FxHashMap::default();
        table.insert(b"AACGTGAT".to_vec(), b"sample_A".to_vec());
        table.insert(b"TGGTGGTA".to_vec(), b"sample_B".to_vec());

        let op = LookupOp::new(
            Label::new(b"seq2.bc1").unwrap(),
            "sample",
            table,
            b"unassigned",
        );

        let reads = vec![make_test_read(b"AACGTGAT")];
        let (result, _) = <LookupOp as GraphNode<NoTrace>>::run_inner(&op, reads).unwrap();
        let reads = result.unwrap();

        let sample = reads[0]
            .data(StrType::Seq(2), InlineString::new(b"bc1"), InlineString::new(b"sample"))
            .unwrap();

        match sample {
            Data::Bytes(v) => assert_eq!(v, b"sample_A"),
            _ => panic!("Expected Bytes data"),
        }
    }

    #[test]
    fn test_lookup_not_found_default() {
        let mut table = FxHashMap::default();
        table.insert(b"AACGTGAT".to_vec(), b"sample_A".to_vec());

        let op = LookupOp::new(
            Label::new(b"seq2.bc1").unwrap(),
            "sample",
            table,
            b"unassigned",
        );

        let reads = vec![make_test_read(b"XXXXXXXX")];
        let (result, _) = <LookupOp as GraphNode<NoTrace>>::run_inner(&op, reads).unwrap();
        let reads = result.unwrap();

        let sample = reads[0]
            .data(StrType::Seq(2), InlineString::new(b"bc1"), InlineString::new(b"sample"))
            .unwrap();

        match sample {
            Data::Bytes(v) => assert_eq!(v, b"unassigned"),
            _ => panic!("Expected Bytes data"),
        }
    }

    #[test]
    fn test_lookup_multiple_reads() {
        let mut table = FxHashMap::default();
        table.insert(b"AACGTGAT".to_vec(), b"sample_A".to_vec());
        table.insert(b"TGGTGGTA".to_vec(), b"sample_B".to_vec());

        let op = LookupOp::new(
            Label::new(b"seq2.bc1").unwrap(),
            "sample",
            table,
            b"unassigned",
        );

        let reads = vec![
            make_test_read(b"AACGTGAT"),
            make_test_read(b"TGGTGGTA"),
            make_test_read(b"UNKNOWN_"),
        ];
        let (result, _) = <LookupOp as GraphNode<NoTrace>>::run_inner(&op, reads).unwrap();
        let reads = result.unwrap();

        let get_sample = |read: &Read| -> Vec<u8> {
            match read
                .data(StrType::Seq(2), InlineString::new(b"bc1"), InlineString::new(b"sample"))
                .unwrap()
            {
                Data::Bytes(v) => v.clone(),
                _ => panic!("Expected Bytes"),
            }
        };

        assert_eq!(get_sample(&reads[0]), b"sample_A");
        assert_eq!(get_sample(&reads[1]), b"sample_B");
        assert_eq!(get_sample(&reads[2]), b"unassigned");
    }
}
