use crate::graph::*;
use crate::inline_string::InlineString;

static COMP_LUT: [u8; 256] = {
    let mut l = [0u8; 256];
    let mut i = 0;

    while i < l.len() {
        l[i] = i as u8;
        i += 1;
    }

    l[b'A' as usize] = b'T';
    l[b'C' as usize] = b'G';
    l[b'G' as usize] = b'C';
    l[b'T' as usize] = b'A';
    l[b'a' as usize] = b't';
    l[b'c' as usize] = b'g';
    l[b'g' as usize] = b'c';
    l[b't' as usize] = b'a';
    l
};

/// Reverse-complement a byte slice in-place.
#[inline]
fn reverse_complement(seq: &mut [u8]) {
    seq.reverse();
    for b in seq.iter_mut() {
        *b = COMP_LUT[*b as usize];
    }
}

pub struct TryOrientationOp<T: Trace = NoTrace> {
    /// The subgraph to execute (geometry matching nodes, no InputOp).
    inner: Graph<T>,
    /// Which read's sequence to reverse-complement on retry (e.g., 1 for seq1).
    read_idx: u8,
    /// Attribute name to store the orientation result (e.g., "ori").
    attr_name: InlineString,
    /// Temporary attribute name for batch-index tracking.
    batch_idx_attr: InlineString,
}

impl<T: Trace> TryOrientationOp<T> {
    const NAME: &'static str = "TryOrientationOp";

    /// Try running the inner graph on each read. If a read is dropped (match
    /// failure), reverse-complement `seq{read_idx}` and retry.
    ///
    /// The orientation that succeeded is stored as `Data::Bytes` (`b"fw"` or
    /// `b"rc"`) in the attribute `attr_name` on the wildcard mapping of
    /// `seq{read_idx}`.
    pub fn new(inner: Graph<T>, read_idx: u8, attr_name: impl AsRef<[u8]>) -> Self {
        Self {
            inner,
            read_idx,
            attr_name: InlineString::new(attr_name.as_ref()),
            batch_idx_attr: InlineString::new(b"_batch_idx"),
        }
    }
}

impl<T: Trace> GraphNode<T> for TryOrientationOp<T> {
    fn run(&self, reads: Option<Vec<Read>>, trace: &T) -> Result<(Option<Vec<Read>>, bool)> {
        let start = trace.start(&reads);
        let Some(reads) = reads else {
            panic!("Expected some reads!")
        };

        let seq_type = StrType::Seq(self.read_idx);
        let wildcard = InlineString::new(b"*");
        let n = reads.len();

        // Step 1: Tag each read with a batch index and clone the batch.
        let mut tagged = reads;
        for (i, read) in tagged.iter_mut().enumerate() {
            *read
                .data_mut(seq_type, wildcard, self.batch_idx_attr)
                .unwrap_or_else(|e| panic!("Error in {}: {e}", Self::NAME)) = Data::Int(i as isize);
        }
        let clones = tagged.clone();

        // Step 2: Forward pass.
        let (fw_result, done) = self.inner.run_one(Some(tagged), trace)?;

        if done {
            let res = (fw_result, done);
            trace.add(self.name(), start, &res.0);
            return Ok(res);
        }

        let fw_survivors = fw_result.unwrap_or_default();

        // Collect surviving batch indices from forward pass.
        let mut fw_survived = vec![false; n];
        for read in &fw_survivors {
            if let Ok(data) = read.data(seq_type, wildcard, self.batch_idx_attr) {
                if let Ok(idx) = data.as_int() {
                    let idx = idx as usize;
                    if idx < n {
                        fw_survived[idx] = true;
                    }
                }
            }
        }

        // Step 3: Build retry batch from clones at dropped indices.
        // RC the target sequence and reverse its quality scores.
        let mut retry_batch: Vec<Read> = Vec::new();
        for (i, mut read) in clones.into_iter().enumerate() {
            if fw_survived[i] {
                continue;
            }

            // Reverse-complement the sequence for seq{read_idx}.
            if let Some(sm) = read.str_mappings_mut(seq_type) {
                let seq = sm.string_mut();
                reverse_complement(seq);

                if let Some(qual) = sm.qual_mut() {
                    qual.reverse();
                }
            }

            retry_batch.push(read);
        }

        // Step 4: RC pass (only if there are reads to retry).
        let rc_survivors = if !retry_batch.is_empty() {
            let (rc_result, _) = self.inner.run_one(Some(retry_batch), trace)?;
            rc_result.unwrap_or_default()
        } else {
            Vec::new()
        };

        // Step 5: Set orientation attributes and remove batch_idx.
        let fw_val = Data::Bytes(b"fw".to_vec());
        let rc_val = Data::Bytes(b"rc".to_vec());

        let mut all_survivors: Vec<(usize, Read)> =
            Vec::with_capacity(fw_survivors.len() + rc_survivors.len());

        for mut read in fw_survivors {
            let idx = read
                .data(seq_type, wildcard, self.batch_idx_attr)
                .expect("_batch_idx must be set on all reads by Step 1")
                .as_int()
                .expect("_batch_idx must be an Int") as usize;

            *read
                .data_mut(seq_type, wildcard, self.attr_name)
                .unwrap_or_else(|e| panic!("Error in {}: {e}", Self::NAME)) = fw_val.clone();

            all_survivors.push((idx, read));
        }

        for mut read in rc_survivors {
            let idx = read
                .data(seq_type, wildcard, self.batch_idx_attr)
                .expect("_batch_idx must be set on all reads by Step 1")
                .as_int()
                .expect("_batch_idx must be an Int") as usize;

            *read
                .data_mut(seq_type, wildcard, self.attr_name)
                .unwrap_or_else(|e| panic!("Error in {}: {e}", Self::NAME)) = rc_val.clone();

            all_survivors.push((idx, read));
        }

        // Step 6: Sort by original batch index to preserve ordering.
        all_survivors.sort_by_key(|(idx, _)| *idx);

        // Step 7: Remove internal _batch_idx attribute from survivors.
        let final_reads: Vec<Read> = all_survivors
            .into_iter()
            .map(|(_, mut r)| {
                r.remove_data(seq_type, wildcard, &self.batch_idx_attr);
                r
            })
            .collect();

        let res = if final_reads.is_empty() {
            (None, false)
        } else {
            (Some(final_reads), false)
        };

        trace.add(self.name(), start, &res.0);
        Ok(res)
    }

    fn required_names(&self) -> &[LabelOrAttr] {
        &[]
    }

    fn name(&self) -> &'static str {
        Self::NAME
    }
}
