use crate::graph::*;

pub struct RetainOp {
    required_names: Vec<LabelOrAttr>,
    selector_expr: Expr,
}

impl RetainOp {
    const NAME: &'static str = "RetainOp";

    /// Retain only the reads where the selector expression evaluates to true and discard the rest.
    pub fn new(selector_expr: impl Into<Expr>) -> Self {
        let selector_expr = selector_expr.into();
        Self {
            required_names: selector_expr.required_names(),
            selector_expr,
        }
    }
}

impl<T: Trace> GraphNode<T> for RetainOp {
    fn run_inner(&self, mut reads: Vec<Read>) -> Result<(Option<Vec<Read>>, bool)> {
        let mut error = None;
        reads.retain(|read| {
            match self.selector_expr.eval_bool(read) {
                Ok(keep) => keep,
                Err(e) => {
                    error = Some(Error::NameError {
                        source: e,
                        read: read.clone(),
                        context: Self::NAME,
                    });
                    false // drop read if error? or stop?
                }
            }
        });
        
        if let Some(e) = error {
            return Err(e);
        }
        
        if reads.is_empty() {
            Ok((None, false))
        } else {
            Ok((Some(reads), false))
        }
    }

    fn required_names(&self) -> &[LabelOrAttr] {
        &self.required_names
    }

    fn name(&self) -> &'static str {
        Self::NAME
    }
}
