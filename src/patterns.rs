use std::borrow::Cow;

use crate::errors::*;
use crate::expr::*;
use crate::inline_string::*;
use crate::read::*;

pub struct Patterns {
    pattern_name: Option<InlineString>,
    multimatch_name: Option<InlineString>,
    attr_names: Vec<InlineString>,
    patterns: Vec<Pattern>,
}

impl Patterns {
    pub fn from_strs(patterns: impl IntoIterator<Item = impl AsRef<[u8]>>) -> Self {
        Self {
            pattern_name: None,
            multimatch_name: None,
            attr_names: Vec::new(),
            patterns: patterns
                .into_iter()
                .map(|v| Pattern::from_literal(v.as_ref(), Vec::new()))
                .collect(),
        }
    }

    pub fn from_exprs(patterns: impl IntoIterator<Item = Expr>) -> Self {
        Self {
            pattern_name: None,
            multimatch_name: None,
            attr_names: Vec::new(),
            patterns: patterns
                .into_iter()
                .map(|v| Pattern::from_expr(v, Vec::new()))
                .collect(),
        }
    }

    pub fn new(
        patterns: impl IntoIterator<Item = Pattern>,
        attr_names: impl IntoIterator<Item = impl AsRef<[u8]>>,
    ) -> Self {
        let attr_names = attr_names
            .into_iter()
            .map(|v| InlineString::new(v.as_ref()))
            .collect::<Vec<_>>();
        let patterns = patterns.into_iter().collect::<Vec<_>>();

        for p in &patterns {
            assert_eq!(
                p.attrs().len(),
                attr_names.len(),
                "Each pattern must have the same number of associated attributes!"
            );
        }

        Self {
            pattern_name: None,
            multimatch_name: None,
            attr_names,
            patterns,
        }
    }

    pub fn with_pattern_name(mut self, pattern_name: impl AsRef<[u8]>) -> Self {
        self.pattern_name = Some(InlineString::new(pattern_name.as_ref()));
        self
    }

    pub fn with_multimatch_name(mut self, multimatch_name: impl AsRef<[u8]>) -> Self {
        self.multimatch_name = Some(InlineString::new(multimatch_name.as_ref()));
        self
    }

    pub fn pattern_name(&self) -> Option<InlineString> {
        self.pattern_name
    }

    pub fn multimatch_name(&self) -> Option<InlineString> {
        self.multimatch_name
    }

    pub fn attr_names(&self) -> &[InlineString] {
        &self.attr_names
    }

    pub fn patterns(&self) -> &[Pattern] {
        &self.patterns
    }

    pub fn iter_literals(&self) -> impl Iterator<Item = (usize, &[u8])> {
        self.patterns.iter().enumerate().filter_map(|(i, p)| {
            if let Pattern::Literal { bytes, .. } = p {
                Some((i, bytes.as_slice()))
            } else {
                None
            }
        })
    }

    pub fn iter_exprs(&self) -> impl Iterator<Item = (usize, &Expr)> {
        self.patterns.iter().enumerate().filter_map(|(i, p)| {
            if let Pattern::Expr { expr, .. } = p {
                Some((i, expr))
            } else {
                None
            }
        })
    }
}

pub enum Pattern {
    Literal { bytes: Vec<u8>, attrs: Vec<Data> },
    Expr { expr: Expr, attrs: Vec<Data> },
}

impl Pattern {
    pub fn from_literal(bytes: &[u8], attrs: Vec<Data>) -> Self {
        Self::Literal {
            bytes: bytes.to_owned(),
            attrs,
        }
    }

    pub fn from_expr(mut expr: Expr, attrs: Vec<Data>) -> Self {
        if expr.optimize() {
            let temp = Read::new();
            let bytes = expr
                .eval_bytes(&temp, false)
                .unwrap_or_else(|e| panic!("{e}"))
                .into_owned();
            Self::Literal { bytes, attrs }
        } else {
            Self::Expr { expr, attrs }
        }
    }

    pub fn get<'a>(&'a self, read: &'a Read) -> std::result::Result<Cow<'a, [u8]>, NameError> {
        use Pattern::*;
        match self {
            Literal { bytes, .. } => Ok(Cow::Borrowed(bytes)),
            Expr { expr, .. } => Ok(expr.eval_bytes(read, false)?),
        }
    }

    pub fn attrs(&self) -> &[Data] {
        use Pattern::*;
        match self {
            Literal { attrs, .. } => attrs,
            Expr { attrs, .. } => attrs,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_patterns_from_strs() {
        let p = Patterns::from_strs(vec!["ACGT", "TGCA"]);
        assert_eq!(p.patterns().len(), 2);
        assert!(p.pattern_name().is_none());
        assert!(p.multimatch_name().is_none());
        assert!(p.attr_names().is_empty());
    }

    #[test]
    fn test_patterns_with_pattern_name() {
        let p = Patterns::from_strs(vec!["ACGT"]).with_pattern_name("match_name");
        assert!(p.pattern_name().is_some());
    }

    #[test]
    fn test_patterns_with_multimatch_name() {
        let p = Patterns::from_strs(vec!["ACGT"]).with_multimatch_name("multi");
        assert!(p.multimatch_name().is_some());
    }

    #[test]
    fn test_patterns_iter_literals() {
        let p = Patterns::from_strs(vec!["ACGT", "TGCA", "AAAA"]);
        let literals: Vec<_> = p.iter_literals().collect();
        assert_eq!(literals.len(), 3);
        assert_eq!(literals[0].0, 0);
        assert_eq!(literals[0].1, b"ACGT");
        assert_eq!(literals[1].1, b"TGCA");
        assert_eq!(literals[2].1, b"AAAA");
    }

    #[test]
    fn test_pattern_from_literal() {
        let p = Pattern::from_literal(b"ACGT", vec![Data::Int(1)]);
        assert_eq!(p.attrs().len(), 1);
        let read = Read::new();
        let bytes = p.get(&read).unwrap();
        assert_eq!(bytes.as_ref(), b"ACGT");
    }

    #[test]
    fn test_pattern_attrs() {
        let p = Pattern::from_literal(b"ACGT", vec![Data::Bool(true), Data::Int(42)]);
        assert_eq!(p.attrs().len(), 2);
    }

    #[test]
    fn test_patterns_new_with_attrs() {
        let patterns = vec![
            Pattern::from_literal(b"ACGT", vec![Data::Int(1)]),
            Pattern::from_literal(b"TGCA", vec![Data::Int(2)]),
        ];
        let p = Patterns::new(patterns, vec!["score"]);
        assert_eq!(p.attr_names().len(), 1);
        assert_eq!(p.patterns().len(), 2);
    }

    #[test]
    fn test_patterns_from_exprs() {
        let expr = Expr::from(b"ACGT".to_vec());
        let p = Patterns::from_exprs(vec![expr]);
        assert_eq!(p.patterns().len(), 1);
        // Constant expression should be optimized to a literal
        let literals: Vec<_> = p.iter_literals().collect();
        assert_eq!(literals.len(), 1);
    }
}
