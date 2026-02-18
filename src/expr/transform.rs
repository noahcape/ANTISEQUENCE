use crate::errors::*;
use crate::expr::{Attr, Label, LabelOrAttr};
use crate::parse_utils::*;

#[derive(Debug, Clone)]
pub struct TransformExpr {
    before: Vec<Label>,
    after: Vec<Option<LabelOrAttr>>,
}

impl TransformExpr {
    pub fn new(
        before: impl IntoIterator<Item = Label>,
        after: impl IntoIterator<Item = Option<LabelOrAttr>>,
    ) -> Self {
        Self {
            before: before.into_iter().collect(),
            after: after.into_iter().collect(),
        }
    }

    /// Parse a byte string to get a transform expression.
    pub fn from_bytes(expr: impl AsRef<[u8]>) -> Result<Self> {
        let (before, after) = parse(expr.as_ref())?;
        Ok(Self { before, after })
    }

    pub fn check_size(&self, before_size: usize, after_size: usize, context: &'static str) {
        assert_eq!(
            before_size,
            self.before.len(),
            "Number of labels before the \"->\" must be {}, only found {} labels for {}",
            before_size,
            self.before.len(),
            context
        );
        assert_eq!(after_size, self.after.len(), "Number of labels or attributes after the \"->\" must be {}, only found {} labels or attributes for {}", after_size, self.after.len(), context);
    }

    pub fn check_same_str_type(&self, context: &'static str) {
        let str_type = self.before[0].str_type;
        assert!(
            self.before.iter().all(|l| l.str_type == str_type),
            "String types before the \"->\" must be the same for {}",
            context
        );
        assert!(
            self.after.iter().all(|label_or_attr| label_or_attr
                .as_ref()
                .map(|l| l.str_type() == str_type)
                .unwrap_or(true)),
            "String types after the \"->\" must be the same for {}",
            context
        );
    }

    pub fn before(&self, i: usize) -> Label {
        self.before[i].clone()
    }

    pub fn after_label(&self, i: usize, context: &'static str) -> Option<Label> {
        self.after[i].clone().map(|a| if let LabelOrAttr::Label(l) = a {
            l
        } else {
            panic!("Expected type.label after the \"->\" in the transform expression for {context}")
        })
    }

    pub fn after_attr(&self, i: usize, context: &'static str) -> Option<Attr> {
        self.after[i].clone().map(|a| if let LabelOrAttr::Attr(a) = a {
            a
        } else {
            panic!("Expected type.label.attr after the \"->\" in the transform expression for {context}")
        })
    }
}

fn parse(expr: &[u8]) -> Result<(Vec<Label>, Vec<Option<LabelOrAttr>>)> {
    let split_idx = expr
        .windows(2)
        .position(|w| w == b"->")
        .ok_or_else(|| Error::Parse {
            string: utf8(expr),
            context: utf8(expr),
            reason: "missing \"->\"",
        })?;
    let before_str = expr[..split_idx].to_owned();
    let after_str = expr[split_idx + 2..].to_owned();

    let before = before_str
        .split(|&b| b == b',')
        .map(Label::new)
        .collect::<Result<Vec<_>>>()?;

    let after = after_str
        .split(|&b| b == b',')
        .map(|s| {
            let s = trim_ascii_whitespace(s).ok_or_else(|| Error::InvalidName {
                string: utf8(s),
                context: utf8(expr),
            })?;

            if s.iter().all(|&c| c == b'_') {
                Ok(None)
            } else {
                Ok(Some(LabelOrAttr::new(s)?))
            }
        })
        .collect::<Result<Vec<_>>>()?;

    Ok((before, after))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_transform_expr_from_bytes() {
        let te = TransformExpr::from_bytes("seq1.a -> seq1.b").unwrap();
        let b = te.before(0);
        assert_eq!(b.label, crate::inline_string::InlineString::new(b"a"));
        let a = te.after_label(0, "test");
        assert!(a.is_some());
    }

    #[test]
    fn test_transform_expr_multiple() {
        let te = TransformExpr::from_bytes("seq1.a, seq1.b -> seq1.c, seq1.d").unwrap();
        te.check_size(2, 2, "test");
    }

    #[test]
    fn test_transform_expr_discard() {
        let te = TransformExpr::from_bytes("seq1.a -> _").unwrap();
        let a = te.after_label(0, "test");
        assert!(a.is_none());
    }

    #[test]
    fn test_transform_expr_with_attr() {
        let te = TransformExpr::from_bytes("seq1.a -> seq1.b.attr").unwrap();
        let a = te.after_attr(0, "test");
        assert!(a.is_some());
    }

    #[test]
    fn test_transform_expr_missing_arrow() {
        let result = TransformExpr::from_bytes("seq1.a seq1.b");
        assert!(result.is_err());
    }

    #[test]
    fn test_check_same_str_type() {
        let te = TransformExpr::from_bytes("seq1.a, seq1.b -> seq1.c, seq1.d").unwrap();
        te.check_same_str_type("test");
    }

    #[test]
    fn test_transform_new() {
        let before = vec![Label::new(b"seq1.a").unwrap()];
        let after = vec![Some(LabelOrAttr::Label(Label::new(b"seq1.b").unwrap()))];
        let te = TransformExpr::new(before, after);
        te.check_size(1, 1, "test");
    }
}
