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
