use std::{
    fs::File,
    io::{BufRead, BufReader},
};

use crate::iter::*;

pub struct FilterReads<R: Reads> {
    reads: R,
    selector_expr: SelectorExpr,
    label: Label,
    attr: Option<Attr>,
    allow_list: Vec<Vec<u8>>,
    mismatch: usize,
}

pub fn parse_allowlist(filename: impl AsRef<Path>) -> Vec<Vec<u8>> {
    let file = File::open(filename.as_ref()).expect("no such file");

    BufReader::new(file)
        .lines()
        .map(|l| {
            let seq = l.expect("Could not parse line");
            seq.as_bytes().to_vec()
        })
        .collect::<Vec<_>>()
}

impl<R: Reads> FilterReads<R> {
    pub fn new(
        reads: R,
        selector_expr: SelectorExpr,
        transform_expr: TransformExpr,
        allow_list: impl AsRef<Path>,
        mismatch: usize,
    ) -> Self {
        transform_expr.check_size(1, 1, "checking length in bounds");
        transform_expr.check_same_str_type("checking length in bounds");

        let allow_list = parse_allowlist(allow_list);

        Self {
            reads,
            selector_expr,
            label: transform_expr.before()[0].clone(),
            attr: transform_expr.after()[0].clone().map(|a| match a {
                LabelOrAttr::Attr(a) => a,
                _ => panic!("Expected type.label.attr after the \"->\" in the transform expression when checking length in bounds"),
            }),
            allow_list,
            mismatch
        }
    }
}

impl<R: Reads> Reads for FilterReads<R> {
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
                read.filter(
                    self.label.str_type,
                    self.label.label,
                    attr.clone(),
                    self.allow_list.clone(),
                    self.mismatch,
                )
                .map_err(|e| Error::NameError {
                    source: e,
                    read: read.clone(),
                    context: "filtering reads",
                })?;
            }
        }

        Ok(reads)
    }

    fn finish(&mut self) -> Result<()> {
        self.reads.finish()
    }
}
