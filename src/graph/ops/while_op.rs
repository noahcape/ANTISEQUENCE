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
    fn run(&self, read: Option<Read>, trace: &T) -> Result<(Option<Read>, bool)> {
        let start = trace.start(&read);
        let Some(mut read) = read else {
            panic!("Expected some read!")
        };

        while self
            .cond_expr
            .eval_bool(&read)
            .map_err(|e| Error::NameError {
                source: e,
                read: read.clone(),
                context: Self::NAME,
            })?
        {
            let (r, done) = self.graph.run_one(Some(read), trace)?;

            if done {
                trace.add(self.name(), start, &r);
                return Ok((r, done));
            }

            if let Some(r) = r {
                read = r;
            } else {
                trace.add(self.name(), start, &r);
                return Ok((r, done));
            }
        }

        let read = Some(read);
        trace.add(self.name(), start, &read);
        Ok((read, false))
    }

    fn required_names(&self) -> &[LabelOrAttr] {
        &self.required_names
    }

    fn name(&self) -> &'static str {
        Self::NAME
    }

    fn finish(&self) -> Result<bool> {
        Ok(true)
    }
}
