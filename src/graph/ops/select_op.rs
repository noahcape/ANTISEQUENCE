use crate::graph::*;

pub struct SelectOp {
    required_names: Vec<LabelOrAttr>,
    selector_expr: Expr,
    graph: Graph,
}

impl SelectOp {
    const NAME: &'static str = "SelectOp";

    /// Run the graph only on reads where the selector expression evaluates to true.
    pub fn new(selector_expr: impl Into<Expr>, graph: Graph) -> Self {
        let selector_expr = selector_expr.into();
        let required_names = selector_expr.required_names();
        Self {
            required_names,
            selector_expr,
            graph,
        }
    }
}

impl GraphNode for SelectOp {
    fn run(&self, read: Option<Read>) -> Result<(Option<Read>, bool)> {
        let Some(read) = read else {
            panic!("Expected some read!")
        };

        if self
            .selector_expr
            .eval_bool(&read)
            .map_err(|e| Error::NameError {
                source: e,
                read: read.clone(),
                context: Self::NAME,
            })?
        {
            self.graph.run_one(Some(read))
        } else {
            Ok((Some(read), false))
        }
    }

    fn required_names(&self) -> &[LabelOrAttr] {
        &self.required_names
    }

    fn name(&self) -> &'static str {
        Self::NAME
    }
}
