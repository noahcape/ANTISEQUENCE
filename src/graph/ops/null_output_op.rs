use crate::graph::*;

pub struct NullOutputOp;

impl Default for NullOutputOp {
    fn default() -> Self {
        Self::new()
    }
}

impl NullOutputOp {
    const NAME: &'static str = "NullOutputOp";

    pub fn new() -> Self {
        Self
    }
}

impl<T: Trace> GraphNode<T> for NullOutputOp {
    fn run_inner(&self, reads: Vec<Read>) -> Result<(Option<Vec<Read>>, bool)> {
        Ok((Some(reads), false))
    }

    fn required_names(&self) -> &[LabelOrAttr] {
        &[]
    }

    fn name(&self) -> &'static str {
        Self::NAME
    }
}
