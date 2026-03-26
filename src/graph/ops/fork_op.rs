use crate::graph::*;

pub struct ForkOp<T: Trace = NoTrace> {
    graph: Graph<T>,
}

impl<T: Trace> ForkOp<T> {
    const NAME: &'static str = "ForkOp";

    /// Clone each read and run the clone through the specified graph, while leaving
    /// the original read unchanged.
    pub fn new(graph: Graph<T>) -> Self {
        Self { graph }
    }
}

impl<T: Trace> GraphNode<T> for ForkOp<T> {
    fn run(&self, reads: Option<Vec<Read>>, trace: &T) -> Result<(Option<Vec<Read>>, bool)> {
        let start = trace.start(&reads);
        let Some(reads) = reads else {
            panic!("Expected some reads!")
        };
        self.graph.run_one(Some(reads.clone()), trace)?;
        let reads = Some(reads);
        trace.add(self.name(), start, &reads);
        Ok((reads, false))
    }

    fn required_names(&self) -> &[LabelOrAttr] {
        &[]
    }

    fn name(&self) -> &'static str {
        Self::NAME
    }
}
