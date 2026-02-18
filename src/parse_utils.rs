use crate::errors::*;
use crate::expr::*;

pub fn trim_ascii_whitespace(b: &[u8]) -> Option<&[u8]> {
    let start = b.iter().position(|&c| !c.is_ascii_whitespace())?;
    let end = b.iter().rposition(|&c| !c.is_ascii_whitespace())?;
    Some(&b[start..=end])
}

pub fn check_valid_name(b: &[u8]) -> Option<&[u8]> {
    for &c in b {
        match c {
            b'a'..=b'z' | b'A'..=b'Z' | b'0'..=b'9' | b'_' | b'*' => (),
            _ => return None,
        }
    }

    Some(b)
}

pub fn parse_fmt_expr(expr: &[u8]) -> Result<Vec<Expr>> {
    let mut res = Vec::new();
    let mut curr = Vec::new();
    let mut escape = false;
    let mut in_label = false;

    for &c in expr {
        match c {
            b'{' if !escape => {
                if in_label {
                    Err(Error::Parse {
                        string: utf8(expr),
                        context: utf8(expr),
                        reason: "cannot have nested braces",
                    })?;
                }
                res.push(Expr::from(curr.clone()));
                in_label = true;
                curr.clear();
            }
            b'}' if !escape => {
                if !in_label {
                    Err(Error::Parse {
                        string: utf8(expr),
                        context: utf8(expr),
                        reason: "unbalanced braces",
                    })?;
                }

                let label = trim_ascii_whitespace(&curr).ok_or_else(|| Error::InvalidName {
                    string: utf8(&curr),
                    context: utf8(expr),
                })?;
                let e = match LabelOrAttr::new(label)? {
                    LabelOrAttr::Label(label) => Expr::from(label),
                    LabelOrAttr::Attr(attr) => Expr::from(attr),
                };
                res.push(e);
                in_label = false;
                curr.clear();
            }
            b'\\' if !escape => escape = true,
            _ => {
                escape = false;
                curr.push(c);
            }
        }
    }

    if !curr.is_empty() {
        if in_label {
            Err(Error::Parse {
                string: utf8(expr),
                context: utf8(expr),
                reason: "unbalanced braces",
            })?;
        }
        res.push(Expr::from(curr));
    }

    Ok(res)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_trim_ascii_whitespace() {
        assert_eq!(
            trim_ascii_whitespace(b"  hello  "),
            Some(b"hello".as_slice())
        );
        assert_eq!(trim_ascii_whitespace(b"hello"), Some(b"hello".as_slice()));
        assert_eq!(trim_ascii_whitespace(b"  "), None);
        assert_eq!(trim_ascii_whitespace(b""), None);
        assert_eq!(trim_ascii_whitespace(b" a "), Some(b"a".as_slice()));
    }

    #[test]
    fn test_check_valid_name() {
        assert!(check_valid_name(b"hello").is_some());
        assert!(check_valid_name(b"test_123").is_some());
        assert!(check_valid_name(b"*").is_some());
        assert!(check_valid_name(b"ABC").is_some());
        assert!(check_valid_name(b"has space").is_none());
        assert!(check_valid_name(b"has.dot").is_none());
        assert!(check_valid_name(b"has-dash").is_none());
    }

    #[test]
    fn test_parse_fmt_expr_plain_text() {
        let result = parse_fmt_expr(b"hello world").unwrap();
        assert_eq!(result.len(), 1);
    }

    #[test]
    fn test_parse_fmt_expr_with_label() {
        let result = parse_fmt_expr(b"prefix{seq1.label}suffix").unwrap();
        assert_eq!(result.len(), 3);
    }

    #[test]
    fn test_parse_fmt_expr_nested_braces_error() {
        let result = parse_fmt_expr(b"{{nested}}");
        assert!(result.is_err());
    }

    #[test]
    fn test_parse_fmt_expr_unbalanced_brace() {
        let result = parse_fmt_expr(b"unbalanced}");
        assert!(result.is_err());
    }

    #[test]
    fn test_parse_fmt_expr_unclosed_brace() {
        let result = parse_fmt_expr(b"{unclosed");
        assert!(result.is_err());
    }

    #[test]
    fn test_parse_fmt_expr_empty() {
        let result = parse_fmt_expr(b"").unwrap();
        assert!(result.is_empty());
    }

    #[test]
    fn test_parse_fmt_expr_escape() {
        let result = parse_fmt_expr(b"hello\\{world").unwrap();
        assert_eq!(result.len(), 1);
    }
}
