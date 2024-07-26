use crate::graph::*;

pub struct ForkOp {
    graph: Graph,
}

impl ForkOp {
    const NAME: &'static str = "ForkOp";

    /// Clone each read and run the clone through the specified graph, while leaving
    /// the original read unchanged.
    pub fn new(graph: Graph) -> Self {
        Self { graph }
    }
}

impl GraphNode for ForkOp {
    fn run(&self, read: Option<Read>) -> Result<(Option<Read>, bool)> {
        let Some(read) = read else {
            panic!("Expected some read!")
        };
        self.graph.run_one(Some(read.clone()))?;
        Ok((Some(read), false))
    }

    fn required_names(&self) -> &[LabelOrAttr] {
        &[]
    }

    fn name(&self) -> &'static str {
        Self::NAME
    }
}
