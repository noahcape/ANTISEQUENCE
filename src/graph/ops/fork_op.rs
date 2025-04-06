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
    fn run(&self, read: Option<Read>, trace: &T) -> Result<(Option<Read>, bool)> {
        let start = trace.start(&read);
        let Some(read) = read else {
            panic!("Expected some read!")
        };
        self.graph.run_one(Some(read.clone()), trace)?;
        let read = Some(read);
        trace.add(self.name(), start, &read);
        Ok((read, false))
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
