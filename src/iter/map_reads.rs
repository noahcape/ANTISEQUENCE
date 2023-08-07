use serde::Deserialize;
use std::{collections::HashMap, path::Path};

use crate::iter::*;

pub struct MapReads<R: Reads> {
    reads: R,
    selector_expr: SelectorExpr,
    label: Label,
    attr: Option<Attr>,
    seq_map: HashMap<Vec<u8>, Vec<u8>>,
    mismatch: usize,
}

#[derive(Debug, Deserialize)]
struct BCMapRecord {
    oligo_dt: String,
    rand_hex: String,
}

pub fn generate_map<'a>(seq_map: impl AsRef<Path>) -> HashMap<Vec<u8>, Vec<u8>> {
    let mut rdr = csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .comment(Some(b'#'))
        .has_headers(false)
        .from_path(seq_map.as_ref())
        .expect(format!("Could not open file {:?}", seq_map.as_ref()).as_str()); // create a custom error for this

    let mut hm: HashMap<Vec<u8>, Vec<u8>> = HashMap::new();

    for result in rdr.deserialize() {
        let record: BCMapRecord = result.expect("could not deseralize map record");

        hm.insert(
            record.rand_hex.as_bytes().to_vec(),
            record.oligo_dt.as_bytes().to_vec(),
        );
    }

    hm
}

impl<R: Reads> MapReads<R> {
    pub fn new(
        reads: R,
        selector_expr: SelectorExpr,
        transform_expr: TransformExpr,
        seq_map: impl AsRef<Path>,
        mismatch: usize,
    ) -> Self {
        transform_expr.check_size(1, 1, "checking length in bounds");
        transform_expr.check_same_str_type("checking length in bounds");

        let seq_map = generate_map(seq_map);

        Self {
            reads,
            selector_expr,
            label: transform_expr.before()[0].clone(),
            attr: transform_expr.after()[0].clone().map(|a| match a {
                LabelOrAttr::Attr(a) => a,
                _ => panic!("Expected type.label.attr after the \"->\" in the transform expression when checking length in bounds"),
            }),
            seq_map,
            mismatch
        }
    }
}

impl<R: Reads> Reads for MapReads<R> {
    fn next_chunk(&self) -> Result<Vec<Read>> {
        let mut reads = self.reads.next_chunk()?;

        for read in reads.iter_mut() {
            if !(self
                .selector_expr
                .matches(read)
                .map_err(|e| Error::NameError {
                    source: e,
                    read: read.clone(),
                    context: "mapping reads",
                })?)
            {
                continue;
            }

            if let Some(attr) = &self.attr {
                read.map(
                    self.label.str_type,
                    self.label.label,
                    attr.clone(),
                    self.seq_map.clone(),
                    self.mismatch,
                )
                .map_err(|e| Error::NameError {
                    source: e,
                    read: read.clone(),
                    context: "mapping read",
                })?;
            }
        }

        Ok(reads)
    }

    fn finish(&mut self) -> Result<()> {
        self.reads.finish()
    }
}
