use crate::graph::*;

pub struct SelectOp<T: Trace = NoTrace> {
    required_names: Vec<LabelOrAttr>,
    selector_expr: Expr,
    graph: Graph<T>,
}

impl<T: Trace> SelectOp<T> {
    const NAME: &'static str = "SelectOp";

    /// Run the graph only on reads where the selector expression evaluates to true.
    pub fn new(selector_expr: impl Into<Expr>, graph: Graph<T>) -> Self {
        let selector_expr = selector_expr.into();
        let required_names = selector_expr.required_names();
        Self {
            required_names,
            selector_expr,
            graph,
        }
    }
}

impl<T: Trace> GraphNode<T> for SelectOp<T> {
    fn run(&self, reads: Option<Vec<Read>>, trace: &T) -> Result<(Option<Vec<Read>>, bool)> {
        let start = trace.start(&reads);
        let Some(reads) = reads else {
            panic!("Expected some reads!")
        };

        // Preserve original ordering: record which indices pass the selector,
        // batch those through the inner graph once, then merge back in order.
        let total = reads.len();
        let mut passing_positions = Vec::new();
        let mut passing_reads = Vec::new();
        let mut failing_reads = Vec::new();

        for (idx, read) in reads.into_iter().enumerate() {
            if self
                .selector_expr
                .eval_bool(&read)
                .map_err(|e| Error::NameError {
                    source: e,
                    read: read.clone(),
                    context: Self::NAME,
                })? {
                passing_positions.push(idx);
                passing_reads.push(read);
            } else {
                failing_reads.push(read);
            }
        }

        // Run the inner graph only on the selected subset (if any).
        let (processed_passing, done) = if !passing_reads.is_empty() {
            let (r, d) = self.graph.run_one(Some(passing_reads), trace)?;
            (r.unwrap_or_default(), d)
        } else {
            (Vec::new(), false)
        };

        // Merge processed passing reads and untouched failing reads back into
        // a single Vec<Read> in the original order.
        let mut res_reads = Vec::with_capacity(total);
        let mut pass_iter = processed_passing.into_iter();
        let mut fail_iter = failing_reads.into_iter();
        let mut pass_pos_idx = 0usize;

        for i in 0..total {
            if pass_pos_idx < passing_positions.len() && passing_positions[pass_pos_idx] == i {
                if let Some(r) = pass_iter.next() {
                    res_reads.push(r);
                }
                pass_pos_idx += 1;
            } else if let Some(r) = fail_iter.next() {
                res_reads.push(r);
            }
        }

        let final_res = if res_reads.is_empty() { None } else { Some(res_reads) };

        trace.add(self.name(), start, &final_res);
        Ok((final_res, done))
    }

    fn required_names(&self) -> &[LabelOrAttr] {
        &self.required_names
    }

    fn name(&self) -> &'static str {
        Self::NAME
    }
}
