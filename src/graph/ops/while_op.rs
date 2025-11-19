use crate::graph::*;

pub struct WhileOp<T: Trace = NoTrace> {
    required_names: Vec<LabelOrAttr>,
    cond_expr: Expr,
    graph: Graph<T>,
}

impl<T: Trace> WhileOp<T> {
    const NAME: &'static str = "WhileOp";

    /// Run a read through the graph multiple times, while the condition expression evaluates to true.
    pub fn new(cond_expr: impl Into<Expr>, graph: Graph<T>) -> Self {
        let cond_expr = cond_expr.into();
        let required_names = cond_expr.required_names();
        Self {
            required_names,
            cond_expr,
            graph,
        }
    }
}

impl<T: Trace> GraphNode<T> for WhileOp<T> {
    fn run(&self, reads: Option<Vec<Read>>, trace: &T) -> Result<(Option<Vec<Read>>, bool)> {
        let start = trace.start(&reads);
        let Some(mut current_batch) = reads else {
            panic!("Expected some reads!")
        };

        let mut final_results = Vec::with_capacity(current_batch.len());
        let mut done_global = false;

        // Loop until no reads are left to process
        while !current_batch.is_empty() {
            let mut passing_reads = Vec::with_capacity(current_batch.len());
            let mut failing_reads = Vec::with_capacity(current_batch.len());

            for read in current_batch {
                if self.cond_expr.eval_bool(&read).map_err(|e| Error::NameError {
                    source: e,
                    read: read.clone(),
                    context: Self::NAME,
                })? {
                    passing_reads.push(read);
                } else {
                    failing_reads.push(read);
                }
            }

            // Failing reads are done
            final_results.extend(failing_reads);

            if passing_reads.is_empty() {
                break;
            }

            // Run passing reads
            let (res_opt, done) = self.graph.run_one(Some(passing_reads), trace)?;
            if done {
                done_global = true;
            }

            if let Some(next_batch) = res_opt {
                current_batch = next_batch;
            } else {
                current_batch = Vec::new();
            }
        }

        let res = if final_results.is_empty() { None } else { Some(final_results) };
        trace.add(self.name(), start, &res);
        Ok((res, done_global))
    }

    fn required_names(&self) -> &[LabelOrAttr] {
        &self.required_names
    }

    fn name(&self) -> &'static str {
        Self::NAME
    }
}
