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
    fn run_inner(&self, read: Read) -> Result<(Option<Read>, bool)> {
        let first_idx = read.first_idx();

        if self.bounds.contains(&first_idx) {
            Ok((Some(read), false))
        } else {
            use std::ops::Bound::*;
            let done = match self.bounds.end_bound() {
                Included(hi) => first_idx > *hi,
                Excluded(hi) => first_idx >= *hi,
                Unbounded => false,
            };
            Ok((None, done))
        }
    }

    fn required_names(&self) -> &[LabelOrAttr] {
        &[]
    }

    fn name(&self) -> &'static str {
        Self::NAME
    }
}
