use std::borrow::Cow;

use crate::errors::*;
use crate::expr::*;
use crate::inline_string::*;
use crate::read::*;

pub struct Patterns {
    pattern_name: Option<InlineString>,
    attr_names: Vec<InlineString>,
    patterns: Vec<Pattern>,
}

impl Patterns {
    pub fn from_strs(patterns: impl IntoIterator<Item = impl AsRef<[u8]>>) -> Self {
        Self {
            pattern_name: None,
            attr_names: Vec::new(),
            patterns: patterns
                .into_iter()
                .map(|v| Pattern::from_literal(v, Vec::new()))
                .collect(),
        }
    }

    pub fn from_exprs(patterns: impl IntoIterator<Item = Expr>) -> Self {
        Self {
            pattern_name: None,
            attr_names: Vec::new(),
            patterns: patterns
                .into_iter()
                .map(|v| Pattern::from_expr(v, Vec::new()))
                .collect(),
        }
    }

    pub fn new(
        pattern_name: impl AsRef<[u8]>,
        attr_names: impl IntoIterator<Item = impl AsRef<[u8]>>,
        patterns: impl IntoIterator<Item = Pattern>,
    ) -> Self {
        Self {
            pattern_name: Some(InlineString::new(pattern_name.as_ref())),
            attr_names: attr_names
                .into_iter()
                .map(|v| InlineString::new(v.as_ref()))
                .collect(),
            patterns: patterns.into_iter().collect(),
        }
    }

    pub fn pattern_name(&self) -> Option<InlineString> {
        self.pattern_name
    }

    pub fn attr_names(&self) -> &[InlineString] {
        &self.attr_names
    }

    pub fn patterns(&self) -> &[Pattern] {
        &self.patterns
    }

    pub fn all_literals(&self) -> bool {
        for p in &self.patterns {
            if let Pattern::Expr { .. } = p {
                return false;
            }
        }

        true
    }

    pub fn iter_literals(&self) -> impl Iterator<Item = &[u8]> {
        self.patterns.iter().map(|p| {
            if let Pattern::Literal { bytes, .. } = p {
                bytes
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
        Self::Literal { bytes: bytes.to_owned(), attrs }
    }

    pub fn from_expr(mut expr: Expr, attrs: Vec<Data>) -> Self {
        if expr.optimize() {
            let temp = Read::new();
            let bytes = expr.eval_bytes(&temp).unwrap_or_else(|e| panic!("{e}")).into_owned();
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
