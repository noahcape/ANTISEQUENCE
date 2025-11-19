use crate::graph::*;

pub struct TakeOp<B: RangeBounds<usize> + Send + Sync> {
    bounds: B,
}

impl<B: RangeBounds<usize> + Send + Sync> TakeOp<B> {
    const NAME: &'static str = "TakeOp";

    /// Take only the reads that have a record index inside the specified bounds.
    pub fn new(bounds: B) -> Self {
        Self { bounds }
    }
}

impl<B: RangeBounds<usize> + Send + Sync, T: Trace> GraphNode<T> for TakeOp<B> {
    fn run_inner(&self, mut reads: Vec<Read>) -> Result<(Option<Vec<Read>>, bool)> {
        let last_idx_in_batch = reads.last().map(|r| r.first_idx());
        
        reads.retain(|read| self.bounds.contains(&read.first_idx()));
        
        use std::ops::Bound::*;
        let done = if let Some(last_idx) = last_idx_in_batch {
             match self.bounds.end_bound() {
                Included(hi) => last_idx > *hi,
                Excluded(hi) => last_idx >= *hi,
                Unbounded => false,
            }
        } else {
            false
        };

        if reads.is_empty() {
            Ok((None, done))
        } else {
            Ok((Some(reads), done))
        }
    }

    fn required_names(&self) -> &[LabelOrAttr] {
        &[]
    }

    fn name(&self) -> &'static str {
        Self::NAME
    }
}
