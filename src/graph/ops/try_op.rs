use crate::graph::*;

pub struct TryOp<T: Trace = NoTrace> {
    try_graph: Graph<T>,
    catch_graph: Graph<T>,
}

impl<T: Trace> TryOp<T> {
    const NAME: &'static str = "TryOp";

    /// Run reads through the try graph, remove the ones that have skipped an operation,
    /// and then run the skipped reads through the catch graph.
    ///
    /// An operation is skipped only if the read does not have a name (label or attribute)
    /// that is required by the operation.
    /// This is useful for specifying a chain of operations where each operation depends on the
    /// labels or attributes produced by the previous operation.
    pub fn new(try_graph: Graph<T>, catch_graph: Graph<T>) -> Self {
        Self {
            try_graph,
            catch_graph,
        }
    }
}

impl<T: Trace> GraphNode<T> for TryOp<T> {
    fn run(&self, read: Option<Read>, trace: &T) -> Result<(Option<Read>, bool)> {
        let start = trace.start(&read);
        let (read, failed, done) = self.try_graph.try_run_one(read, trace)?;

        let res = if !done && read.is_some() && failed {
            let (_, done) = self.catch_graph.run_one(read, trace)?;
            (None, done)
        } else {
            (read, done)
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

    fn finish(&self) -> Result<bool> {
        Ok(true)
    }
}
