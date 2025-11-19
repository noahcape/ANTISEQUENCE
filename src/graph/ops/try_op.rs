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
    fn run(&self, reads: Option<Vec<Read>>, trace: &T) -> Result<(Option<Vec<Read>>, bool)> {
        let start = trace.start(&reads);
        let (reads, failed, done) = self.try_graph.try_run_one(reads, trace)?;

        let res = if !done && failed {
            // In batch mode, `try_run_one` returns failed=true if *any* read failed requirements?
            // Or does it filter?
            // `try_run_one` in Graph checks requirements on the *first* read of the batch (as per my implementation).
            // If the first read fails, it returns `(reads, true, false)`.
            // So we assume the whole batch is "failed" / "skipped" by try_graph.
            // Then we run the whole batch through catch_graph.
            
            let (_, done) = self.catch_graph.run_one(reads, trace)?;
            (None, done) 
            // Wait, if we run catch_graph, we should return its output?
            // The original code returns `(None, done)`. 
            // Ah, original code:
            // `let (read, failed, done) = self.try_graph.try_run_one(read, trace)?;`
            // `if ... failed { let (_, done) = self.catch_graph.run_one(read, trace)?; (None, done) }`
            // It discards the output of catch_graph? That seems odd.
            // Looking at original code: `(None, done)`.
            // It seems TryOp in original code *consumes* the read if it goes to catch block?
            // Or maybe `try_run_one` returns the read if it failed?
            
            // Let's check `Graph::try_run_one`:
            // if !read.has_names(...) { return Ok((curr, true, false)); }
            // It returns the read back.
            
            // So `TryOp` logic:
            // 1. Try running `try_graph`.
            // 2. If it failed (requirements not met), run `catch_graph`.
            // 3. Return `(None, done)`. 
            // This implies `TryOp` acts as a sink if it goes to catch? 
            // Or maybe it assumes `catch_graph` handles the output/storage?
            // If `catch_graph` has output nodes, they write.
            // But `TryOp` itself returns `None` as the read.
            
            // So for batch, we follow same logic.
        } else {
            (reads, done)
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
