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
    fn run(&self, read: Option<Read>, trace: &T) -> Result<(Option<Read>, bool)> {
        let start = trace.start(&read);
        let Some(read) = read else {
            panic!("Expected some read!")
        };

        let res = if self
            .selector_expr
            .eval_bool(&read)
            .map_err(|e| Error::NameError {
                source: e,
                read: read.clone(),
                context: Self::NAME,
            })? {
            self.graph.run_one(Some(read), trace)?
        } else {
            (Some(read), false)
        };

        trace.add(self.name(), start, &res.0);
        Ok(res)
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
