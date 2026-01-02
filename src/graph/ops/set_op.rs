use crate::graph::*;

/// Operation to set a label or attribute on a Read.
///
/// This operation evaluates an expression for each read and assigns the result to
/// a specified label (modifying the sequence/quality) or attribute (metadata).
pub struct SetOp {
    required_names: Vec<LabelOrAttr>,
    label_or_attr: LabelOrAttr,
    expr: Expr,
}

impl SetOp {
    const NAME: &'static str = "SetOp";

    /// Set a labeled interval or attribute to the result of an expression.
    ///
    /// The expression must return a byte string if a labeled interval is being set.
    ///
    /// To generate the quality scores when setting intervals that have corresponding quality
    /// scores, references to intervals in the expression are directly substituted with the
    /// corresponding quality scores of the intervals. For references to byte strings without
    /// quality scores, a sequence of `I`s is used as the quality scores in the expression.
    /// *This naive substitution may lead to unexpected results for complex expressions!*
    ///
    /// If a label is set, then its interval and all other intersecting intervals will be adjusted accordingly
    /// for any shortening or lengthening.
    pub fn new(label_or_attr: impl Into<LabelOrAttr>, expr: impl Into<Expr>) -> Self {
        let label_or_attr = label_or_attr.into();
        let mut expr = expr.into();
        expr.optimize();
        let mut required_names = expr.required_names();
        match &label_or_attr {
            // if label_or_attr is a label, we need to add it to required_names to ensure the label is present
            LabelOrAttr::Label(_) => required_names.push(label_or_attr.clone()),
            // otherwise if label_or_attr is an attribute, we need to add the label it depends on to required_names
            LabelOrAttr::Attr(a) => required_names.push(LabelOrAttr::Label(Label {
                str_type: a.str_type,
                label: a.label,
            })),
        }

        Self {
            required_names,
            label_or_attr,
            expr,
        }
    }
}

impl<T: Trace> GraphNode<T> for SetOp {
    /// Executes the Set operation on a batch of reads.
    ///
    /// This method iterates through the provided reads and modifies them based on the
    /// configured target (`label_or_attr`) and expression.
    ///
    /// # Logic
    /// * **If targeting a Label (Sequence/Quality):**
    ///   1. Evaluates the expression to generate the new sequence bytes.
    ///   2. If the read has quality scores, it re-evaluates the expression in "quality mode"
    ///      to generate corresponding quality scores.
    ///   3. Updates the read's sequence (and quality) at the specified label.
    ///
    /// * **If targeting an Attribute (Metadata):**
    ///   1. Evaluates the expression to compute the new attribute value.
    ///   2. Updates the read's metadata storage with this new value.
    fn run_inner(&self, mut reads: Vec<Read>) -> Result<(Option<Vec<Read>>, bool)> {
        for read in &mut reads {
            match &self.label_or_attr {
                // Case 1: label_or_attr is a label
                LabelOrAttr::Label(label) => {
                    // Evaluate the expression to get the new byte sequence
                    let new_bytes = self
                        .expr
                        .eval_bytes(&read, false)
                        .map_err(|e| Error::NameError {
                            source: e,
                            read: read.clone(),
                            context: Self::NAME,
                        })?
                        .into_owned();

                    let str_mappings =
                        read.str_mappings(label.str_type)
                            .ok_or_else(|| Error::NameError {
                                source: NameError::NotInRead(Name::StrType(label.str_type)),
                                read: read.clone(),
                                context: Self::NAME,
                            })?;

                    // If quality scores exist, evaluate the expression for quality scores as well
                    // This generates corresponding quality scores for the new sequence
                    // (e.g., if you concatenate two sequences, it concatenates their quality scores).
                    if str_mappings.qual().is_some() {
                        let new_qual = self
                            .expr
                            .eval_bytes(&read, true)
                            .map_err(|e| Error::NameError {
                                source: e,
                                read: read.clone(),
                                context: Self::NAME,
                            })?
                            .into_owned();

                        // Update the read content (sequence and quality)
                        read.set(label.str_type, label.label, &new_bytes, Some(&new_qual))
                            .map_err(|e| Error::NameError {
                                source: e,
                                read: read.clone(),
                                context: Self::NAME,
                            })?;
                    } else {
                        // Update only the sequence content
                        read.set(label.str_type, label.label, &new_bytes, None)
                            .map_err(|e| Error::NameError {
                                source: e,
                                read: read.clone(),
                                context: Self::NAME,
                            })?;
                    }
                }
                LabelOrAttr::Attr(attr) => {
                    // Evaluate expression for attribute value
                    let new_val = self.expr.eval(&read, false).map_err(|e| Error::NameError {
                        source: e,
                        read: read.clone(),
                        context: Self::NAME,
                    })?;

                    // Update the attribute in the read's metadata
                    // panic to make borrow checker happy
                    *read
                        .data_mut(attr.str_type, attr.label, attr.attr)
                        .unwrap_or_else(|e| panic!("Error in {}: {e}", Self::NAME)) = new_val.into();
                }
            }
        }

        Ok((Some(reads), false))
    }

    fn required_names(&self) -> &[LabelOrAttr] {
        &self.required_names
    }

    fn name(&self) -> &'static str {
        Self::NAME
    }
}
