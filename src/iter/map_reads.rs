use needletail::bitkmer::BitNuclKmer;
use serde::Deserialize;
use std::collections::HashMap;

use crate::iter::*;

pub struct MapReads<R: Reads> {
    reads: R,
    selector_expr: SelectorExpr,
    label: Label,
    attr: Option<Attr>,
    seq_map: String,
    mismatch: usize,
}

#[derive(Debug, Deserialize)]
struct BCMapRecord {
    oligo_dt: String,
    rand_hex: String,
}

pub fn generate_maps(
    seq_map: String,
    k: usize,
) -> (HashMap<u64, String>, HashMap<u64, Vec<u64>>) {
    let mut rdr = csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .comment(Some(b'#'))
        .has_headers(false)
        .from_path(seq_map.clone())
        .expect(format!("Could not open file {}", seq_map).as_str()); // create a custom error for this

    let mut hm = HashMap::new();
    let mut kmer_hm = HashMap::new();

    for result in rdr.deserialize() {
        let record: BCMapRecord = result.expect("could not deseralize map record");

        let mut kmers: Vec<u64> = Vec::new();

        if let Some((_, (rh, _), _)) = BitNuclKmer::new(
            record.rand_hex.as_bytes(),
            record.rand_hex.len() as u8,
            false,
        )
        .next()
        {
            for i in 0..record.rand_hex.len() - k {
                if let Some((_, (kmer, _), _)) =
                    BitNuclKmer::new(record.rand_hex[i..i + k].as_bytes(), k as u8, false).next()
                {
                    kmers.push(kmer)
                }
            }

            kmer_hm.insert(rh, kmers);
            hm.insert(rh, record.oligo_dt.clone());
        }
    }

    (hm, kmer_hm)
}

pub fn edit_distance(target: u64, t_len: usize, query: u64, q_len: usize) -> usize {
    if t_len == 0 {
        return q_len;
    }

    if q_len == 0 {
        return t_len;
    }

    if target & 3 == query & 3 {
        edit_distance(target >> 2, t_len - 1, query >> 2, q_len - 1)
    } else {
        1 + [
            edit_distance(target >> 2, t_len - 1, query, q_len),
            edit_distance(target, t_len, query >> 2, q_len - 1),
            edit_distance(target >> 2, t_len - 1, query >> 2, q_len - 1),
        ]
        .into_iter()
        .min()
        .unwrap()
    }
}

impl<R: Reads> MapReads<R> {
    pub fn new(
        reads: R,
        selector_expr: SelectorExpr,
        transform_expr: TransformExpr,
        seq_map: String,
        mismatch: usize,
    ) -> Self {
        transform_expr.check_size(1, 1, "checking length in bounds");
        transform_expr.check_same_str_type("checking length in bounds");

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
