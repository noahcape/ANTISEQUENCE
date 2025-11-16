use crate::graph::*;

pub struct NullOutputOp;

impl NullOutputOp {
    const NAME: &'static str = "NullOutputOp";

    pub fn new() -> Self {
        Self
    }
}

impl<T: Trace> GraphNode<T> for NullOutputOp {
    fn run_inner(&self, read: Read) -> Result<(Option<Read>, bool)> {
        Ok((Some(read), false))
    }

    fn required_names(&self) -> &[LabelOrAttr] {
        &[]
    }

    fn name(&self) -> &'static str {
        Self::NAME
    }
}
