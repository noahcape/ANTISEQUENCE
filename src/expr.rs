mod transform;
pub use transform::*;

mod expr_node;
pub use expr_node::*;

use crate::errors::*;
use crate::inline_string::*;
use crate::parse_utils::*;
use crate::read::*;

#[derive(Debug, Clone, PartialEq, Hash)]
pub struct Label {
    pub str_type: StrType,
    pub label: InlineString,
}

#[derive(Debug, Clone, PartialEq, Hash)]
pub struct Attr {
    pub str_type: StrType,
    pub label: InlineString,
    pub attr: InlineString,
}

#[derive(Debug, Clone, PartialEq)]
pub enum LabelOrAttr {
    Label(Label),
    Attr(Attr),
}

impl Label {
    pub fn new(s: &[u8]) -> Result<Self> {
        let split = s.split(|&b| b == b'.').collect::<Vec<_>>();

        match split.as_slice() {
            &[str_type, label] => {
                let str_type =
                    trim_ascii_whitespace(str_type).ok_or_else(|| Error::InvalidName {
                        string: utf8(str_type),
                        context: utf8(s),
                    })?;
                let label = trim_ascii_whitespace(label).ok_or_else(|| Error::InvalidName {
                    string: utf8(label),
                    context: utf8(s),
                })?;
                let label = check_valid_name(label).ok_or_else(|| Error::InvalidName {
                    string: utf8(label),
                    context: utf8(s),
                })?;

                Ok(Self {
                    str_type: StrType::new(str_type)?,
                    label: InlineString::new(label),
                })
            }
            _ => Err(Error::Parse {
                string: utf8(s),
                context: utf8(s),
                reason: "expected type.label",
            }),
        }
    }
}

impl Attr {
    pub fn new(s: &[u8]) -> Result<Self> {
        let split = s.split(|&b| b == b'.').collect::<Vec<_>>();

        match split.as_slice() {
            &[str_type, label, attr] => {
                let str_type =
                    trim_ascii_whitespace(str_type).ok_or_else(|| Error::InvalidName {
                        string: utf8(str_type),
                        context: utf8(s),
                    })?;
                let label = trim_ascii_whitespace(label).ok_or_else(|| Error::InvalidName {
                    string: utf8(label),
                    context: utf8(s),
                })?;
                let label = check_valid_name(label).ok_or_else(|| Error::InvalidName {
                    string: utf8(label),
                    context: utf8(s),
                })?;

                let attr = trim_ascii_whitespace(attr).ok_or_else(|| Error::InvalidName {
                    string: utf8(attr),
                    context: utf8(s),
                })?;
                let attr = check_valid_name(attr).ok_or_else(|| Error::InvalidName {
                    string: utf8(attr),
                    context: utf8(s),
                })?;

                Ok(Self {
                    str_type: StrType::new(str_type)?,
                    label: InlineString::new(label),
                    attr: InlineString::new(attr),
                })
            }
            _ => Err(Error::Parse {
                string: utf8(s),
                context: utf8(s),
                reason: "expected type.label.attr",
            }),
        }
    }
}

impl LabelOrAttr {
    pub fn new(s: &[u8]) -> Result<Self> {
        let count = s.iter().filter(|&&c| c == b'.').count();

        match count {
            1 => Ok(LabelOrAttr::Label(Label::new(s)?)),
            2 => Ok(LabelOrAttr::Attr(Attr::new(s)?)),
            _ => Err(Error::Parse {
                string: utf8(s),
                context: utf8(s),
                reason: "expected type.label or type.label.attr",
            }),
        }
    }

    pub fn str_type(&self) -> StrType {
        match self {
            LabelOrAttr::Label(l) => l.str_type,
            LabelOrAttr::Attr(a) => a.str_type,
        }
    }

    pub fn label(&self) -> InlineString {
        match self {
            LabelOrAttr::Label(l) => l.label,
            LabelOrAttr::Attr(a) => a.label,
        }
    }
}

impl From<Label> for LabelOrAttr {
    fn from(label: Label) -> Self {
        LabelOrAttr::Label(label)
    }
}

impl From<Attr> for LabelOrAttr {
    fn from(attr: Attr) -> Self {
        LabelOrAttr::Attr(attr)
    }
}

/// Create a transform expression.
#[macro_export]
macro_rules! tr {
    ($($t:tt)+) => {
        {
            let s = stringify!($($t)+);
            $crate::expr::TransformExpr::from_bytes(s.as_bytes())
                .unwrap_or_else(|e| panic!("Error constructing transform expression:\n{e}\non line {} column {} in file {}", line!(), column!(), file!()))
        }
    };
}

/// Create a label by parsing a byte string of the form `type.label`.
pub fn label(s: impl AsRef<[u8]>) -> Label {
    Label::new(s.as_ref()).unwrap_or_else(|e| panic!("Error creating label:\n{e}"))
}

/// Create an attribute by parsing a byte string of the form `type.label.attr`.
pub fn attr(s: impl AsRef<[u8]>) -> Attr {
    Attr::new(s.as_ref()).unwrap_or_else(|e| panic!("Error creating attr:\n{e}"))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_label_new_valid() {
        let l = Label::new(b"seq1.test").unwrap();
        assert_eq!(l.str_type, StrType::Seq(1));
        assert_eq!(l.label, InlineString::new(b"test"));
    }

    #[test]
    fn test_label_new_star() {
        let l = Label::new(b"seq1.*").unwrap();
        assert_eq!(l.label, InlineString::new(b"*"));
    }

    #[test]
    fn test_label_new_name_type() {
        let l = Label::new(b"name1.readname").unwrap();
        assert_eq!(l.str_type, StrType::Name(1));
    }

    #[test]
    fn test_label_new_invalid_format() {
        assert!(Label::new(b"notype").is_err());
        assert!(Label::new(b"too.many.parts").is_err());
    }

    #[test]
    fn test_label_new_invalid_str_type() {
        assert!(Label::new(b"invalid.label").is_err());
    }

    #[test]
    fn test_label_new_with_whitespace() {
        let l = Label::new(b" seq1 . test ").unwrap();
        assert_eq!(l.str_type, StrType::Seq(1));
        assert_eq!(l.label, InlineString::new(b"test"));
    }

    #[test]
    fn test_attr_new_valid() {
        let a = Attr::new(b"seq1.label.attr").unwrap();
        assert_eq!(a.str_type, StrType::Seq(1));
        assert_eq!(a.label, InlineString::new(b"label"));
        assert_eq!(a.attr, InlineString::new(b"attr"));
    }

    #[test]
    fn test_attr_new_invalid() {
        assert!(Attr::new(b"seq1.label").is_err());
        assert!(Attr::new(b"notype").is_err());
    }

    #[test]
    fn test_label_or_attr_label() {
        let la = LabelOrAttr::new(b"seq1.test").unwrap();
        assert!(matches!(la, LabelOrAttr::Label(_)));
        assert_eq!(la.str_type(), StrType::Seq(1));
        assert_eq!(la.label(), InlineString::new(b"test"));
    }

    #[test]
    fn test_label_or_attr_attr() {
        let la = LabelOrAttr::new(b"seq1.label.attr").unwrap();
        assert!(matches!(la, LabelOrAttr::Attr(_)));
        assert_eq!(la.str_type(), StrType::Seq(1));
        assert_eq!(la.label(), InlineString::new(b"label"));
    }

    #[test]
    fn test_label_or_attr_invalid() {
        assert!(LabelOrAttr::new(b"notype").is_err());
        assert!(LabelOrAttr::new(b"a.b.c.d").is_err());
    }

    #[test]
    fn test_label_or_attr_from_label() {
        let l = Label::new(b"seq1.test").unwrap();
        let la: LabelOrAttr = l.into();
        assert!(matches!(la, LabelOrAttr::Label(_)));
    }

    #[test]
    fn test_label_or_attr_from_attr() {
        let a = Attr::new(b"seq1.label.attr").unwrap();
        let la: LabelOrAttr = a.into();
        assert!(matches!(la, LabelOrAttr::Attr(_)));
    }

    #[test]
    fn test_label_helper() {
        let l = label("seq1.test");
        assert_eq!(l.str_type, StrType::Seq(1));
    }

    #[test]
    fn test_attr_helper() {
        let a = attr("seq1.label.attr");
        assert_eq!(a.str_type, StrType::Seq(1));
    }
}
