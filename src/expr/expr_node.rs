use std::borrow::Cow;
use std::marker::{Send, Sync};
use std::ops::{Bound, RangeBounds};

use crate::errors::NameError;
use crate::expr::*;
use crate::read::*;

const UNKNOWN_QUAL: u8 = b'I';

// Default DNA
pub const NUC_MAP: [u8; 4] = [b'A', b'C', b'T', b'G'];

/// One node in an expression tree.
pub struct Expr {
    node: Box<dyn ExprNode + Send + Sync>,
}

macro_rules! binary_fn {
    ($fn_name:ident, $struct_name:ident) => {
        pub fn $fn_name(self, o: impl Into<Expr>) -> Expr {
            Expr {
                node: Box::new($struct_name {
                    left: self,
                    right: o.into(),
                }),
            }
        }
    };
}

macro_rules! unary_fn {
    ($fn_name:ident, $struct_name:ident, $field_name:ident) => {
        pub fn $fn_name(self) -> Expr {
            Expr {
                node: Box::new($struct_name { $field_name: self }),
            }
        }
    };
}

pub fn log4_roundup(n: usize) -> usize {
    std::ops::Div::div((usize::BITS - n.leading_zeros()) as f64, 2.0).ceil() as usize
}

impl Expr {
    binary_fn!(and, AndNode);
    binary_fn!(or, OrNode);
    binary_fn!(xor, XorNode);

    binary_fn!(add, AddNode);
    binary_fn!(sub, SubNode);
    binary_fn!(mul, MulNode);
    binary_fn!(div, DivNode);

    binary_fn!(gt, GtNode);
    binary_fn!(lt, LtNode);
    binary_fn!(ge, GeNode);
    binary_fn!(le, LeNode);
    binary_fn!(eq, EqNode);

    binary_fn!(concat, ConcatNode);

    unary_fn!(not, NotNode, boolean);
    unary_fn!(len, LenNode, string);
    unary_fn!(rev, RevNode, string);

    unary_fn!(int, IntNode, convert);
    unary_fn!(float, FloatNode, convert);
    unary_fn!(bytes, BytesNode, convert);

    pub fn repeat(self, times: impl Into<Expr>) -> Expr {
        Expr {
            node: Box::new(RepeatNode {
                string: self,
                times: times.into(),
            }),
        }
    }

    pub fn in_bounds<E, R: RangeBounds<E>>(
        self,
        range: impl RangeInto<E, R, (Bound<Expr>, Bound<Expr>)>,
    ) -> Expr {
        Expr {
            node: Box::new(InBoundsNode {
                num: self,
                range: range.range_into(),
            }),
        }
    }

    pub fn slice<E, R: RangeBounds<E>>(
        self,
        range: impl RangeInto<E, R, (Bound<Expr>, Bound<Expr>)>,
    ) -> Expr {
        Expr {
            node: Box::new(SliceNode {
                string: self,
                range: range.range_into(),
            }),
        }
    }

    pub fn revcomp_rna(self) -> Expr {
        Expr {
            node: Box::new(RevCompNode {
                string: self,
                is_rna: true,
            }),
        }
    }

    pub fn revcomp(self) -> Expr {
        Expr {
            node: Box::new(RevCompNode {
                string: self,
                is_rna: false,
            }),
        }
    }

    pub fn pad(self, pad_char: Expr, num: Expr, end: End) -> Expr {
        Expr {
            node: Box::new(PadNode {
                string: self,
                pad_char,
                num,
                end,
            }),
        }
    }

    pub fn normalize<E, R: RangeBounds<E>>(
        self,
        range: impl RangeInto<E, R, (Bound<Expr>, Bound<Expr>)>,
    ) -> Expr {
        Expr {
            node: Box::new(NormalizeNode {
                string: self,
                range: range.range_into(),
            }),
        }
    }

    pub fn eval_bool<'a>(&'a self, read: &'a Read) -> std::result::Result<bool, NameError> {
        let res = self.eval(read, false)?;

        if let EvalData::Bool(b) = res {
            Ok(b)
        } else {
            Err(NameError::Type("bool", vec![res.to_data()]))
        }
    }

    pub fn eval_bytes<'a>(
        &'a self,
        read: &'a Read,
        use_qual: bool,
    ) -> std::result::Result<Cow<'a, [u8]>, NameError> {
        let res = self.eval(read, use_qual)?;

        if let EvalData::Bytes(b) = res {
            Ok(b)
        } else {
            Err(NameError::Type("bytes", vec![res.to_data()]))
        }
    }

    pub fn eval<'a>(
        &'a self,
        read: &'a Read,
        use_qual: bool,
    ) -> std::result::Result<EvalData<'a>, NameError> {
        self.node.eval(read, use_qual)
    }

    pub fn required_names(&self) -> Vec<LabelOrAttr> {
        self.node.required_names()
    }
}

pub trait ExprNode {
    fn eval<'a>(
        &'a self,
        read: &'a Read,
        use_qual: bool,
    ) -> std::result::Result<EvalData<'a>, NameError>;
    fn required_names(&self) -> Vec<LabelOrAttr>;
}

macro_rules! bool_binary_ops {
    ($struct_name:ident, $bool_expr:expr) => {
        struct $struct_name {
            left: Expr,
            right: Expr,
        }

        impl ExprNode for $struct_name {
            fn eval<'a>(
                &'a self,
                read: &'a Read,
                use_qual: bool,
            ) -> std::result::Result<EvalData<'a>, NameError> {
                let left = self.left.eval(read, use_qual)?;
                let right = self.right.eval(read, use_qual)?;

                use EvalData::*;
                match (left, right) {
                    (Bool(l), Bool(r)) => Ok($bool_expr(l, r)),
                    (l, r) => Err(NameError::Type("bool", vec![l.to_data(), r.to_data()])),
                }
            }

            fn required_names(&self) -> Vec<LabelOrAttr> {
                let mut res = self.left.required_names();
                res.append(&mut self.right.required_names());
                res
            }
        }
    };
}

bool_binary_ops!(AndNode, |l, r| Bool(l & r));
bool_binary_ops!(OrNode, |l, r| Bool(l | r));
bool_binary_ops!(XorNode, |l, r| Bool(l ^ r));

macro_rules! num_binary_ops {
    ($struct_name:ident, $int_expr:expr, $float_expr:expr) => {
        struct $struct_name {
            left: Expr,
            right: Expr,
        }

        impl ExprNode for $struct_name {
            fn eval<'a>(
                &'a self,
                read: &'a Read,
                use_qual: bool,
            ) -> std::result::Result<EvalData<'a>, NameError> {
                let left = self.left.eval(read, use_qual)?;
                let right = self.right.eval(read, use_qual)?;

                use EvalData::*;
                match (left, right) {
                    (Int(l), Int(r)) => Ok($int_expr(l, r)),
                    (Float(l), Float(r)) => Ok($float_expr(l, r)),
                    (l, r) => Err(NameError::Type(
                        "both int or both float",
                        vec![l.to_data(), r.to_data()],
                    )),
                }
            }

            fn required_names(&self) -> Vec<LabelOrAttr> {
                let mut res = self.left.required_names();
                res.append(&mut self.right.required_names());
                res
            }
        }
    };
}

num_binary_ops!(AddNode, |l, r| Int(l + r), |l, r| Float(l + r));
num_binary_ops!(SubNode, |l, r| Int(l - r), |l, r| Float(l - r));
num_binary_ops!(MulNode, |l, r| Int(l * r), |l, r| Float(l * r));
num_binary_ops!(DivNode, |l, r| Int(l / r), |l, r| Float(l / r));

num_binary_ops!(GtNode, |l, r| Bool(l > r), |l, r| Bool(l > r));
num_binary_ops!(LtNode, |l, r| Bool(l < r), |l, r| Bool(l < r));
num_binary_ops!(GeNode, |l, r| Bool(l >= r), |l, r| Bool(l >= r));
num_binary_ops!(LeNode, |l, r| Bool(l <= r), |l, r| Bool(l <= r));

struct NotNode {
    boolean: Expr,
}

impl ExprNode for NotNode {
    fn eval<'a>(
        &'a self,
        read: &'a Read,
        use_qual: bool,
    ) -> std::result::Result<EvalData<'a>, NameError> {
        let boolean = self.boolean.eval(read, use_qual)?;
        use EvalData::*;
        match boolean {
            Bool(b) => Ok(Bool(!b)),
            b => Err(NameError::Type("bool", vec![b.to_data()])),
        }
    }

    fn required_names(&self) -> Vec<LabelOrAttr> {
        self.boolean.required_names()
    }
}

struct EqNode {
    left: Expr,
    right: Expr,
}

impl ExprNode for EqNode {
    fn eval<'a>(
        &'a self,
        read: &'a Read,
        use_qual: bool,
    ) -> std::result::Result<EvalData<'a>, NameError> {
        let left = self.left.eval(read, use_qual)?;
        let right = self.right.eval(read, use_qual)?;
        use EvalData::*;
        match (left, right) {
            (Int(l), Int(r)) => Ok(Bool(l == r)),
            (Float(l), Float(r)) => Ok(Bool(l == r)),
            (Bool(l), Bool(r)) => Ok(Bool(l == r)),
            (Bytes(l), Bytes(r)) => Ok(Bool(l == r)),
            (l, r) => Err(NameError::Type(
                "both are the same type",
                vec![l.to_data(), r.to_data()],
            )),
        }
    }

    fn required_names(&self) -> Vec<LabelOrAttr> {
        let mut res = self.left.required_names();
        res.append(&mut self.right.required_names());
        res
    }
}

struct NormalizeNode {
    string: Expr,
    range: (Bound<Expr>, Bound<Expr>),
}

impl ExprNode for NormalizeNode {
    fn eval<'a>(
        &'a self,
        read: &'a Read,
        use_qual: bool,
    ) -> std::result::Result<EvalData<'a>, NameError> {
        let string = self.string.eval(read, use_qual)?;

        use EvalData::*;
        let string = match string {
            Bytes(b) => b,
            b => return Err(NameError::Type("bytes", vec![b.to_data()])),
        };

        let mut start_add1 = false;
        let start = match &self.range.0 {
            Bound::Included(s) => s.eval(read, use_qual)?,
            Bound::Excluded(s) => {
                start_add1 = true;
                s.eval(read, use_qual)?
            }
            Bound::Unbounded => {
                return Err(NameError::Type(
                    "inclusive or exclusive bounds",
                    vec![Data::Int(std::isize::MIN)],
                ))
            }
        };

        let mut end_sub1 = false;
        let end = match &self.range.1 {
            Bound::Included(e) => e.eval(read, use_qual)?,
            Bound::Excluded(e) => {
                end_sub1 = true;
                e.eval(read, use_qual)?
            }
            Bound::Unbounded => {
                return Err(NameError::Type(
                    "inclusive or exclusive bounds",
                    vec![Data::Int(std::isize::MAX)],
                ))
            }
        };

        let (start, end) = match (start, end) {
            (Int(mut s), Int(mut e)) => {
                if start_add1 {
                    s += 1;
                }

                if end_sub1 {
                    e -= 1;
                }

                (s, e)
            }
            (s, e) => return Err(NameError::Type("both int", vec![s.to_data(), e.to_data()])),
        };

        if end < string.len() as isize {
            return Err(NameError::Type(
                "range end to exceed length of interval",
                vec![Data::Int(end)],
            ));
        }

        let mut length_diff = end as usize - string.len();
        let extra_len = log4_roundup((end - start + 1) as usize);

        let mut buff = [b'A'].repeat(length_diff);
        let mut variable_seg = [b'0'].repeat(extra_len);

        for i in 0..extra_len {
            let nuc = NUC_MAP.get(length_diff & (usize::MAX & 3)).unwrap();
            length_diff >>= 2;

            variable_seg[i] = *nuc;
        }

        let mut normalized = string.to_vec();
        normalized.append(&mut buff);
        normalized.append(&mut variable_seg);

        Ok(Bytes(Cow::Owned(normalized)))
    }

    fn required_names(&self) -> Vec<LabelOrAttr> {
        let mut res = self.string.required_names();
        match self.range.start_bound() {
            Bound::Included(s) | Bound::Excluded(s) => res.append(&mut s.required_names()),
            _ => (),
        }
        match self.range.end_bound() {
            Bound::Included(e) | Bound::Excluded(e) => res.append(&mut e.required_names()),
            _ => (),
        }
        res
    }
}

struct PadNode {
    string: Expr,
    pad_char: Expr,
    num: Expr,
    end: End,
}

impl ExprNode for PadNode {
    fn eval<'a>(
        &'a self,
        read: &'a Read,
        use_qual: bool,
    ) -> std::result::Result<EvalData<'a>, NameError> {
        let string = self.string.eval(read, use_qual)?;
        let pad_char = self.pad_char.eval(read, use_qual)?;
        let num = self.num.eval(read, use_qual)?;

        use EvalData::*;
        let string = match string {
            Bytes(b) => b,
            b => return Err(NameError::Type("bytes", vec![b.to_data()])),
        };

        let num = match num {
            Int(n) => {
                if string.len() as isize > n {
                    return Err(NameError::Type(
                        "usize longer than interval length",
                        vec![num.to_data()],
                    ));
                }
                n as usize
            }
            n => return Err(NameError::Type("int", vec![n.to_data()])),
        };

        let pad_char = match pad_char {
            Bytes(b) => b,
            b => return Err(NameError::Type("bytes", vec![b.to_data()])),
        };

        let padded = match self.end {
            Left => {
                let mut padded = pad_char.repeat(num - string.len());
                padded.append(&mut string.to_vec());
                padded
            }
            Right => {
                let mut padded = string.to_vec();
                padded.append(&mut pad_char.repeat(num - string.len()));
                padded
            }
        };
        Ok(Bytes(Cow::Owned(padded)))
    }

    fn required_names(&self) -> Vec<LabelOrAttr> {
        self.string.required_names()
    }
}

struct RevCompNode {
    string: Expr,
    is_rna: bool,
}

impl ExprNode for RevCompNode {
    fn eval<'a>(
        &'a self,
        read: &'a Read,
        use_qual: bool,
    ) -> std::result::Result<EvalData<'a>, NameError> {
        let string = self.string.eval(read, use_qual)?;

        use EvalData::*;
        match string {
            Bytes(b) => {
                let b = if self.is_rna {
                    b.into_iter()
                        .map(|e| match e {
                            b'A' => b'U',
                            b'U' => b'A',
                            b'G' => b'C',
                            b'C' => b'G',
                            b => *b,
                        })
                        .collect::<Vec<_>>()
                } else {
                    b.into_iter()
                        .map(|e| match e {
                            b'A' => b'T',
                            b'T' => b'A',
                            b'G' => b'C',
                            b'C' => b'G',
                            b => *b,
                        })
                        .collect::<Vec<_>>()
                };

                Ok(Bytes(Cow::Owned(b)))
            }
            b => Err(NameError::Type("bytes", vec![b.to_data()])),
        }
    }

    fn required_names(&self) -> Vec<LabelOrAttr> {
        self.string.required_names()
    }
}

struct RevNode {
    string: Expr,
}

impl ExprNode for RevNode {
    fn eval<'a>(
        &'a self,
        read: &'a Read,
        use_qual: bool,
    ) -> std::result::Result<EvalData<'a>, NameError> {
        let string = self.string.eval(read, use_qual)?;

        use EvalData::*;
        match string {
            Bytes(b) => Ok(Bytes(Cow::Owned(
                b.into_iter().rev().map(|e| *e).collect::<Vec<_>>(),
            ))),
            b => Err(NameError::Type("bytes", vec![b.to_data()])),
        }
    }

    fn required_names(&self) -> Vec<LabelOrAttr> {
        self.string.required_names()
    }
}

struct LenNode {
    string: Expr,
}

impl ExprNode for LenNode {
    fn eval<'a>(
        &'a self,
        read: &'a Read,
        use_qual: bool,
    ) -> std::result::Result<EvalData<'a>, NameError> {
        let string = self.string.eval(read, use_qual)?;
        use EvalData::*;
        match string {
            Bytes(s) => Ok(Int(s.len() as isize)),
            s => Err(NameError::Type("bytes", vec![s.to_data()])),
        }
    }

    fn required_names(&self) -> Vec<LabelOrAttr> {
        self.string.required_names()
    }
}

struct IntNode {
    convert: Expr,
}

impl ExprNode for IntNode {
    fn eval<'a>(
        &'a self,
        read: &'a Read,
        use_qual: bool,
    ) -> std::result::Result<EvalData<'a>, NameError> {
        let convert = self.convert.eval(read, use_qual)?;
        use EvalData::*;
        match convert {
            Bool(c) => Ok(Int(if c { 1 } else { 0 })),
            Int(c) => Ok(Int(c)),
            Float(c) => Ok(Int(c as isize)),
            Bytes(c) => Ok(Int(std::str::from_utf8(&c)
                .unwrap()
                .parse::<isize>()
                .unwrap_or_else(|e| panic!("{e}")))),
        }
    }

    fn required_names(&self) -> Vec<LabelOrAttr> {
        self.convert.required_names()
    }
}

struct FloatNode {
    convert: Expr,
}

impl ExprNode for FloatNode {
    fn eval<'a>(
        &'a self,
        read: &'a Read,
        use_qual: bool,
    ) -> std::result::Result<EvalData<'a>, NameError> {
        let convert = self.convert.eval(read, use_qual)?;
        use EvalData::*;
        match convert {
            Bool(c) => Ok(Float(if c { 1.0 } else { 0.0 })),
            Int(c) => Ok(Float(c as f64)),
            Float(c) => Ok(Float(c)),
            Bytes(c) => Ok(Float(
                std::str::from_utf8(&c)
                    .unwrap()
                    .parse::<f64>()
                    .unwrap_or_else(|e| panic!("{e}")),
            )),
        }
    }

    fn required_names(&self) -> Vec<LabelOrAttr> {
        self.convert.required_names()
    }
}

struct BytesNode {
    convert: Expr,
}

impl ExprNode for BytesNode {
    fn eval<'a>(
        &'a self,
        read: &'a Read,
        use_qual: bool,
    ) -> std::result::Result<EvalData<'a>, NameError> {
        let convert = self.convert.eval(read, use_qual)?;
        use EvalData::*;
        match convert {
            Bool(c) => Ok(Bytes(Cow::Borrowed(if c { b"true" } else { b"false" }))),
            Int(c) => Ok(Bytes(Cow::Owned(c.to_string().into_bytes()))),
            Float(c) => Ok(Bytes(Cow::Owned(c.to_string().into_bytes()))),
            Bytes(c) => Ok(Bytes(c)),
        }
    }

    fn required_names(&self) -> Vec<LabelOrAttr> {
        self.convert.required_names()
    }
}

struct RepeatNode {
    string: Expr,
    times: Expr,
}

impl ExprNode for RepeatNode {
    fn eval<'a>(
        &'a self,
        read: &'a Read,
        use_qual: bool,
    ) -> std::result::Result<EvalData<'a>, NameError> {
        let string = self.string.eval(read, use_qual)?;
        let times = self.times.eval(read, use_qual)?;
        use EvalData::*;
        match (string, times) {
            (Bytes(s), Int(t)) => Ok(Bytes(Cow::Owned(s.repeat(t as usize)))),
            (s, t) => Err(NameError::Type(
                "bytes and int",
                vec![s.to_data(), t.to_data()],
            )),
        }
    }

    fn required_names(&self) -> Vec<LabelOrAttr> {
        let mut res = self.string.required_names();
        res.append(&mut self.times.required_names());
        res
    }
}

struct ConcatNode {
    left: Expr,
    right: Expr,
}

impl ExprNode for ConcatNode {
    fn eval<'a>(
        &'a self,
        read: &'a Read,
        use_qual: bool,
    ) -> std::result::Result<EvalData<'a>, NameError> {
        let left = self.left.eval(read, use_qual)?;
        let right = self.right.eval(read, use_qual)?;
        use EvalData::*;
        match (left, right) {
            (Bytes(mut l), Bytes(r)) => {
                l.to_mut().extend_from_slice(&r);
                Ok(Bytes(l))
            }
            (l, r) => Err(NameError::Type(
                "both bytes",
                vec![l.to_data(), r.to_data()],
            )),
        }
    }

    fn required_names(&self) -> Vec<LabelOrAttr> {
        let mut res = self.left.required_names();
        res.append(&mut self.right.required_names());
        res
    }
}

pub fn concat_all(nodes: impl IntoIterator<Item = impl Into<Expr>>) -> Expr {
    Expr {
        node: Box::new(ConcatAllNode {
            nodes: nodes.into_iter().map(|e| e.into()).collect(),
        }),
    }
}

struct ConcatAllNode {
    nodes: Vec<Expr>,
}

impl ExprNode for ConcatAllNode {
    fn eval<'a>(
        &'a self,
        read: &'a Read,
        use_qual: bool,
    ) -> std::result::Result<EvalData<'a>, NameError> {
        let mut res = Vec::new();

        for node in &self.nodes {
            let b = node.eval(read, use_qual)?;

            if let EvalData::Bytes(b) = b {
                res.extend_from_slice(&b);
            } else {
                return Err(NameError::Type("bytes", vec![b.to_data()]));
            }
        }

        Ok(EvalData::Bytes(Cow::Owned(res)))
    }

    fn required_names(&self) -> Vec<LabelOrAttr> {
        self.nodes.iter().flat_map(|n| n.required_names()).collect()
    }
}

struct SliceNode {
    string: Expr,
    range: (Bound<Expr>, Bound<Expr>),
}

impl ExprNode for SliceNode {
    fn eval<'a>(
        &'a self,
        read: &'a Read,
        use_qual: bool,
    ) -> std::result::Result<EvalData<'a>, NameError> {
        let string = self.string.eval(read, use_qual)?;

        use EvalData::*;
        match string {
            Bytes(b) => {
                let mut start_add1 = false;
                let start = match &self.range.0 {
                    Bound::Included(s) => s.eval(read, use_qual)?,
                    Bound::Excluded(s) => {
                        start_add1 = true;
                        s.eval(read, use_qual)?
                    }
                    Bound::Unbounded => Int(0),
                };

                let mut end_sub1 = false;
                let end = match &self.range.1 {
                    Bound::Included(e) => e.eval(read, use_qual)?,
                    Bound::Excluded(e) => {
                        end_sub1 = true;
                        e.eval(read, use_qual)?
                    }
                    Bound::Unbounded => Int(b.len() as isize),
                };

                match (start, end) {
                    (Int(mut s), Int(mut e)) => {
                        let len = b.len() as isize;
                        if e < 0 {
                            e = len + e;
                        }

                        if end_sub1 {
                            e -= 1
                        }

                        if start_add1 {
                            s += 1
                        }

                        if s >= 0 && s <= len && e <= len && e >= 0 && s <= e {
                            Ok(EvalData::Bytes(Cow::Owned(
                                b.get(s as usize..e as usize).unwrap().into(),
                            )))
                        } else {
                            Err(NameError::Type(
                                "indices in bound",
                                vec![Data::Int(s), Data::Int(e)],
                            ))
                        }
                    }
                    (s, e) => Err(NameError::Type("all int", vec![s.to_data(), e.to_data()])),
                }
            }
            b => Err(NameError::Type("interval", vec![b.to_data()])),
        }
    }

    fn required_names(&self) -> Vec<LabelOrAttr> {
        let mut res = self.string.required_names();
        match self.range.start_bound() {
            Bound::Included(s) | Bound::Excluded(s) => res.append(&mut s.required_names()),
            _ => (),
        }
        match self.range.end_bound() {
            Bound::Included(e) | Bound::Excluded(e) => res.append(&mut e.required_names()),
            _ => (),
        }
        res
    }
}

struct InBoundsNode {
    num: Expr,
    range: (Bound<Expr>, Bound<Expr>),
}

impl ExprNode for InBoundsNode {
    fn eval<'a>(
        &'a self,
        read: &'a Read,
        use_qual: bool,
    ) -> std::result::Result<EvalData<'a>, NameError> {
        let num = self.num.eval(read, use_qual)?;

        use EvalData::*;
        let mut start_add1 = false;
        let start = match self.range.start_bound() {
            Bound::Included(s) => s.eval(read, use_qual)?,
            Bound::Excluded(s) => {
                start_add1 = true;
                s.eval(read, use_qual)?
            }
            Bound::Unbounded => Int(std::isize::MIN),
        };

        let mut end_sub1 = false;
        let end = match self.range.end_bound() {
            Bound::Included(e) => e.eval(read, use_qual)?,
            Bound::Excluded(e) => {
                end_sub1 = true;
                e.eval(read, use_qual)?
            }
            Bound::Unbounded => Int(std::isize::MAX),
        };

        match (num, start, end) {
            (Int(n), Int(mut s), Int(mut e)) => {
                if start_add1 {
                    s += 1;
                }
                if end_sub1 {
                    e -= 1;
                }
                Ok(Bool(s <= n && n <= e))
            }
            (n, s, e) => Err(NameError::Type(
                "all int",
                vec![n.to_data(), s.to_data(), e.to_data()],
            )),
        }
    }

    fn required_names(&self) -> Vec<LabelOrAttr> {
        let mut res = self.num.required_names();
        match self.range.start_bound() {
            Bound::Included(s) => res.append(&mut s.required_names()),
            Bound::Excluded(s) => res.append(&mut s.required_names()),
            _ => (),
        }
        match self.range.end_bound() {
            Bound::Included(e) => res.append(&mut e.required_names()),
            Bound::Excluded(e) => res.append(&mut e.required_names()),
            _ => (),
        }
        res
    }
}

pub fn fmt_expr(format_str: impl AsRef<[u8]>) -> Expr {
    let exprs = parse_fmt_expr(format_str.as_ref())
        .unwrap_or_else(|e| panic!("Error constructing format expression:\n{e}"));
    concat_all(exprs)
}

impl ExprNode for Label {
    fn eval<'a>(
        &'a self,
        read: &'a Read,
        use_qual: bool,
    ) -> std::result::Result<EvalData<'a>, NameError> {
        Ok(EvalData::Bytes(if use_qual {
            if let Some(qual) = read.substring_qual(self.str_type, self.label)? {
                Cow::Borrowed(qual)
            } else {
                Cow::Owned(vec![
                    UNKNOWN_QUAL;
                    read.mapping(self.str_type, self.label)?.len
                ])
            }
        } else {
            Cow::Borrowed(read.substring(self.str_type, self.label)?)
        }))
    }

    fn required_names(&self) -> Vec<LabelOrAttr> {
        vec![LabelOrAttr::Label(self.clone())]
    }
}

impl ExprNode for Attr {
    fn eval<'a>(
        &'a self,
        read: &'a Read,
        use_qual: bool,
    ) -> std::result::Result<EvalData<'a>, NameError> {
        let d = read.data(self.str_type, self.label, self.attr)?;
        if use_qual {
            if let Data::Bytes(b) = d {
                return Ok(EvalData::Bytes(Cow::Owned(vec![UNKNOWN_QUAL; b.len()])));
            }
        }

        match d {
            Data::Bool(b) => Ok(EvalData::Bool(*b)),
            Data::Int(i) => Ok(EvalData::Int(*i)),
            Data::Float(f) => Ok(EvalData::Float(*f)),
            Data::Bytes(b) => Ok(EvalData::Bytes(Cow::Borrowed(b))),
        }
    }

    fn required_names(&self) -> Vec<LabelOrAttr> {
        vec![LabelOrAttr::Attr(self.clone())]
    }
}

pub fn label_exists(name: impl AsRef<[u8]>) -> Expr {
    Expr {
        node: Box::new(LabelExistsNode {
            label: Label::new(name.as_ref()).unwrap_or_else(|e| panic!("{e}")),
        }),
    }
}

struct LabelExistsNode {
    label: Label,
}

impl ExprNode for LabelExistsNode {
    fn eval<'a>(
        &'a self,
        read: &'a Read,
        _use_qual: bool,
    ) -> std::result::Result<EvalData<'a>, NameError> {
        Ok(EvalData::Bool(
            read.mapping(self.label.str_type, self.label.label).is_ok(),
        ))
    }

    fn required_names(&self) -> Vec<LabelOrAttr> {
        Vec::new()
    }
}

pub fn attr_exists(name: impl AsRef<[u8]>) -> Expr {
    Expr {
        node: Box::new(AttrExistsNode {
            attr: Attr::new(name.as_ref()).unwrap_or_else(|e| panic!("{e}")),
        }),
    }
}

struct AttrExistsNode {
    attr: Attr,
}

impl ExprNode for AttrExistsNode {
    fn eval<'a>(
        &'a self,
        read: &'a Read,
        _use_qual: bool,
    ) -> std::result::Result<EvalData<'a>, NameError> {
        Ok(EvalData::Bool(
            read.data(self.attr.str_type, self.attr.label, self.attr.attr)
                .is_ok(),
        ))
    }

    fn required_names(&self) -> Vec<LabelOrAttr> {
        Vec::new()
    }
}

impl ExprNode for Data {
    fn eval<'a>(
        &'a self,
        _read: &'a Read,
        use_qual: bool,
    ) -> std::result::Result<EvalData<'a>, NameError> {
        if use_qual {
            if let Data::Bytes(b) = self {
                return Ok(EvalData::Bytes(Cow::Owned(vec![UNKNOWN_QUAL; b.len()])));
            }
        }

        match self {
            Data::Bool(b) => Ok(EvalData::Bool(*b)),
            Data::Int(i) => Ok(EvalData::Int(*i)),
            Data::Float(f) => Ok(EvalData::Float(*f)),
            Data::Bytes(b) => Ok(EvalData::Bytes(Cow::Borrowed(b))),
        }
    }

    fn required_names(&self) -> Vec<LabelOrAttr> {
        Vec::new()
    }
}

impl ExprNode for &str {
    fn eval<'a>(
        &'a self,
        _read: &'a Read,
        use_qual: bool,
    ) -> std::result::Result<EvalData<'a>, NameError> {
        if use_qual {
            Ok(EvalData::Bytes(Cow::Owned(vec![UNKNOWN_QUAL; self.len()])))
        } else {
            Ok(EvalData::Bytes(Cow::Borrowed(self.as_bytes())))
        }
    }

    fn required_names(&self) -> Vec<LabelOrAttr> {
        Vec::new()
    }
}

impl ExprNode for String {
    fn eval<'a>(
        &'a self,
        _read: &'a Read,
        use_qual: bool,
    ) -> std::result::Result<EvalData<'a>, NameError> {
        if use_qual {
            Ok(EvalData::Bytes(Cow::Owned(vec![UNKNOWN_QUAL; self.len()])))
        } else {
            Ok(EvalData::Bytes(Cow::Borrowed(self.as_bytes())))
        }
    }

    fn required_names(&self) -> Vec<LabelOrAttr> {
        Vec::new()
    }
}

impl ExprNode for &[u8] {
    fn eval<'a>(
        &'a self,
        _read: &'a Read,
        use_qual: bool,
    ) -> std::result::Result<EvalData<'a>, NameError> {
        if use_qual {
            Ok(EvalData::Bytes(Cow::Owned(vec![UNKNOWN_QUAL; self.len()])))
        } else {
            Ok(EvalData::Bytes(Cow::Borrowed(self)))
        }
    }

    fn required_names(&self) -> Vec<LabelOrAttr> {
        Vec::new()
    }
}

impl ExprNode for Vec<u8> {
    fn eval<'a>(
        &'a self,
        _read: &'a Read,
        use_qual: bool,
    ) -> std::result::Result<EvalData<'a>, NameError> {
        if use_qual {
            Ok(EvalData::Bytes(Cow::Owned(vec![UNKNOWN_QUAL; self.len()])))
        } else {
            Ok(EvalData::Bytes(Cow::Borrowed(&self)))
        }
    }

    fn required_names(&self) -> Vec<LabelOrAttr> {
        Vec::new()
    }
}

macro_rules! impl_expr_node {
    ($type_name:ident, $eval_data_variant:ident) => {
        impl ExprNode for $type_name {
            fn eval<'a>(
                &'a self,
                _read: &'a Read,
                _use_qual: bool,
            ) -> std::result::Result<EvalData<'a>, NameError> {
                Ok(EvalData::$eval_data_variant(*self as _))
            }

            fn required_names(&self) -> Vec<LabelOrAttr> {
                Vec::new()
            }
        }
    };
}

impl_expr_node!(usize, Int);
impl_expr_node!(isize, Int);
impl_expr_node!(u64, Int);
impl_expr_node!(i64, Int);
impl_expr_node!(u32, Int);
impl_expr_node!(i32, Int);
impl_expr_node!(bool, Bool);
impl_expr_node!(f64, Float);
impl_expr_node!(f32, Float);

impl<E: ExprNode + Send + Sync + 'static> From<E> for Expr {
    fn from(v: E) -> Self {
        Expr { node: Box::new(v) }
    }
}

pub trait FromRange<E, R: RangeBounds<E>> {
    fn from_range(r: R) -> Self;
}

pub trait RangeInto<E, R: RangeBounds<E>, T> {
    fn range_into(self) -> T;
}

impl<E, R: RangeBounds<E>, T: FromRange<E, R>> RangeInto<E, R, T> for R {
    fn range_into(self) -> T {
        T::from_range(self)
    }
}

impl<E: ExprNode + Send + Sync + Copy + 'static, R: RangeBounds<E>> FromRange<E, R>
    for (Bound<Expr>, Bound<Expr>)
{
    fn from_range(r: R) -> Self {
        (
            r.start_bound().map(|&b| Expr::from(b)),
            r.end_bound().map(|&b| Expr::from(b)),
        )
    }
}

macro_rules! impl_from_range_expr {
    ($type_name:ty, $r:ident, $tuple:expr) => {
        impl FromRange<Expr, $type_name> for (Bound<Expr>, Bound<Expr>) {
            fn from_range($r: $type_name) -> Self {
                $tuple
            }
        }
    };
}

impl_from_range_expr!(
    std::ops::Range<Expr>,
    r,
    (Bound::Included(r.start), Bound::Excluded(r.end))
);
impl_from_range_expr!(
    std::ops::RangeFrom<Expr>,
    r,
    (Bound::Included(r.start), Bound::Unbounded)
);
impl_from_range_expr!(
    std::ops::RangeFull,
    _r,
    (Bound::Unbounded, Bound::Unbounded)
);
impl_from_range_expr!(std::ops::RangeInclusive<Expr>, r, {
    let (s, e) = r.into_inner();
    (Bound::Included(s), Bound::Included(e))
});
impl_from_range_expr!(
    std::ops::RangeTo<Expr>,
    r,
    (Bound::Unbounded, Bound::Excluded(r.end))
);
impl_from_range_expr!(
    std::ops::RangeToInclusive<Expr>,
    r,
    (Bound::Unbounded, Bound::Included(r.end))
);

#[derive(Debug)]
pub enum EvalData<'a> {
    Bool(bool),
    Int(isize),
    Float(f64),
    Bytes(Cow<'a, [u8]>),
}

impl<'a> EvalData<'a> {
    pub fn to_data(self) -> Data {
        match self {
            EvalData::Bool(b) => Data::Bool(b),
            EvalData::Int(i) => Data::Int(i),
            EvalData::Float(f) => Data::Float(f),
            EvalData::Bytes(b) => Data::Bytes(b.into_owned()),
        }
    }
}
