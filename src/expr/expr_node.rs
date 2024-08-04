use std::borrow::Cow;
use std::marker::{Send, Sync};
use std::ops::{Bound, RangeBounds};

use crate::errors::NameError;
use crate::expr::*;
use crate::read::*;

const UNKNOWN_QUAL: u8 = b'I';

static NUC: [u8; 4] = [b'A', b'C', b'G', b'T'];
static COMP_LUT: [u8; 256] = {
    let mut l = [0u8; 256];
    let mut i = 0;

    while i < l.len() {
        l[i] = i as u8;
        i += 1;
    }

    l[b'A' as usize] = b'T';
    l[b'C' as usize] = b'G';
    l[b'G' as usize] = b'C';
    l[b'T' as usize] = b'A';
    l
};

pub fn log4_roundup(n: usize) -> usize {
    (n.ilog2() + 1).div_ceil(2) as usize
}

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
    unary_fn!(revcomp, RevCompNode, string);

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

    pub fn pad(self, pad_char: Expr, len: Expr, end: End) -> Expr {
        Expr {
            node: Box::new(PadNode {
                string: self,
                pad_char,
                len,
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
        expect_bool(res)
    }

    pub fn eval_int<'a>(&'a self, read: &'a Read) -> std::result::Result<isize, NameError> {
        let res = self.eval(read, false)?;
        expect_int(res)
    }

    pub fn eval_bytes<'a>(
        &'a self,
        read: &'a Read,
        use_qual: bool,
    ) -> std::result::Result<Cow<'a, [u8]>, NameError> {
        let res = self.eval(read, use_qual)?;
        expect_bytes(res)
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

    /// Apply optimization passes.
    ///
    /// Returns whether the result is just a constant.
    pub fn optimize(&mut self) -> bool {
        let constant = self.propagate_const();
        constant
    }

    fn propagate_const(&mut self) -> bool {
        if !self.required_names().is_empty() {
            return false;
        }

        // no read-specific dependencies
        let temp = Read::new();
        let data = self.eval(&temp, false).unwrap_or_else(|e| panic!("{e}"));
        self.node = Box::new(Data::from(data));
        true
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
                    (l, r) => Err(NameError::Type("bool", vec![l.into(), r.into()])),
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
                        vec![l.into(), r.into()],
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
        Ok(EvalData::Bool(!(expect_bool(boolean)?)))
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
                vec![l.into(), r.into()],
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
        let string = expect_bytes(self.string.eval(read, use_qual)?)?;
        let range = map_eval_range(&self.range, read, use_qual, |b| b as usize)?;
        let (start, end) = range_inclusive_exclusive(&range, std::usize::MAX);

        if end == std::usize::MAX {
            return Err(NameError::Other(
                "the end bound must be bounded for normalization",
            ));
        }
        if end < string.len() {
            return Err(NameError::Other(
                "the string length is greater than the end bound",
            ));
        }

        let pad_len = end - string.len();
        let var_len = log4_roundup(end - start);

        let mut normalized = string.into_owned();
        let pad_char = if use_qual { UNKNOWN_QUAL } else { NUC[0] };
        normalized.extend(std::iter::repeat(pad_char).take(pad_len));

        if use_qual {
            normalized.extend(std::iter::repeat(UNKNOWN_QUAL).take(var_len));
        } else {
            normalized.extend(
                (0..var_len).map(|i| unsafe { *NUC.as_ptr().add((pad_len >> (i * 2)) & 0b11) }),
            );
        }

        Ok(EvalData::Bytes(Cow::Owned(normalized)))
    }

    fn required_names(&self) -> Vec<LabelOrAttr> {
        let mut res = self.string.required_names();
        self.range
            .start_bound()
            .map(|s| res.append(&mut s.required_names()));
        self.range
            .end_bound()
            .map(|e| res.append(&mut e.required_names()));
        res
    }
}

struct PadNode {
    string: Expr,
    pad_char: Expr,
    len: Expr,
    end: End,
}

impl ExprNode for PadNode {
    fn eval<'a>(
        &'a self,
        read: &'a Read,
        use_qual: bool,
    ) -> std::result::Result<EvalData<'a>, NameError> {
        let string = expect_bytes(self.string.eval(read, use_qual)?)?;
        let pad_char = expect_bytes(self.pad_char.eval(read, use_qual)?)?;
        let len = expect_int(self.len.eval(read, use_qual)?)? as usize;

        if string.len() > len {
            return Err(NameError::Other(
                "the length of the string is greater than the padded length",
            ));
        }
        if pad_char.len() != 1 {
            return Err(NameError::Other(
                "the length of padding character must be 1",
            ));
        }

        let padded = match self.end {
            Left => {
                let mut padded = pad_char.repeat(len - string.len());
                padded.extend_from_slice(&string);
                Cow::Owned(padded)
            }
            Right => {
                let mut padded = string.to_owned();
                padded
                    .to_mut()
                    .extend(std::iter::repeat(pad_char[0]).take(len - string.len()));
                padded
            }
        };
        Ok(EvalData::Bytes(padded))
    }

    fn required_names(&self) -> Vec<LabelOrAttr> {
        let mut res = self.string.required_names();
        res.append(&mut self.pad_char.required_names());
        res.append(&mut self.len.required_names());
        res
    }
}

struct RevCompNode {
    string: Expr,
}

impl ExprNode for RevCompNode {
    fn eval<'a>(
        &'a self,
        read: &'a Read,
        use_qual: bool,
    ) -> std::result::Result<EvalData<'a>, NameError> {
        let string = expect_bytes(self.string.eval(read, use_qual)?)?;
        string.to_mut().reverse();

        if !use_qual {
            string.to_mut()
                .iter_mut()
                .for_each(|c| { *c = unsafe { *COMP_LUT.as_ptr().add(*c as usize) }; });
        }

        Ok(EvalData::Bytes(string))
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
        let string = expect_bytes(self.string.eval(read, use_qual)?)?;
        string.to_mut().reverse();
        Ok(EvalData::Bytes(string))
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
        Ok(EvalData::Int(expect_bytes(string)?.len() as isize))
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
            (s, t) => Err(NameError::Type("bytes and int", vec![s.into(), t.into()])),
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
            (l, r) => Err(NameError::Type("both bytes", vec![l.into(), r.into()])),
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
            res.extend_from_slice(&expect_bytes(b)?);
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
        let string = expect_bytes(self.string.eval(read, use_qual)?)?;
        let len = string.len() as isize;
        let range = map_eval_range(&self.range, read, use_qual, |b| {
            if b < 0 {
                (len + b) as usize
            } else {
                b as usize
            }
        })?;

        match string {
            Cow::Owned(mut s) => {
                let (start, end) = range_inclusive_exclusive(&range, len as usize);
                s.copy_within(range, 0);
                s.truncate(end - start);
                Ok(EvalData::Bytes(Cow::Owned(s)))
            }
            Cow::Borrowed(s) => Ok(EvalData::Bytes(Cow::Borrowed(&s[range]))),
        }
    }

    fn required_names(&self) -> Vec<LabelOrAttr> {
        let mut res = self.string.required_names();
        self.range
            .start_bound()
            .map(|s| res.append(&mut s.required_names()));
        self.range
            .end_bound()
            .map(|e| res.append(&mut e.required_names()));
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
        let num = expect_int(self.num.eval(read, use_qual)?)?;
        let range = map_eval_range(&self.range, read, use_qual, |b| b)?;
        Ok(EvalData::Bool(range.contains(&num)))
    }

    fn required_names(&self) -> Vec<LabelOrAttr> {
        let mut res = self.num.required_names();
        self.range
            .start_bound()
            .map(|s| res.append(&mut s.required_names()));
        self.range
            .end_bound()
            .map(|e| res.append(&mut e.required_names()));
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

macro_rules! impl_expr_node_slice {
    ($type_name:ty) => {
        impl ExprNode for $type_name {
            fn eval<'a>(
                &'a self,
                _read: &'a Read,
                use_qual: bool,
            ) -> std::result::Result<EvalData<'a>, NameError> {
                let b: &[u8] = self.as_ref();

                if use_qual {
                    Ok(EvalData::Bytes(Cow::Owned(vec![UNKNOWN_QUAL; b.len()])))
                } else {
                    Ok(EvalData::Bytes(Cow::Borrowed(b)))
                }
            }

            fn required_names(&self) -> Vec<LabelOrAttr> {
                Vec::new()
            }
        }
    };
}

macro_rules! impl_expr_node {
    ($type_name:ty, $eval_data_variant:ident) => {
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

impl_expr_node_slice!(&[u8]);
impl_expr_node_slice!(Vec<u8>);
impl_expr_node_slice!(&str);
impl_expr_node_slice!(String);
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

impl<'a> From<EvalData<'a>> for Data {
    fn from(v: EvalData<'a>) -> Self {
        match v {
            EvalData::Bool(b) => Data::Bool(b),
            EvalData::Int(i) => Data::Int(i),
            EvalData::Float(f) => Data::Float(f),
            EvalData::Bytes(b) => Data::Bytes(b.into_owned()),
        }
    }
}

macro_rules! impl_expect_type {
    ($fn_name:ident, $type_name:ty, $eval_data_name:pat, $v:ident, $str:expr) => {
        fn $fn_name<'a>(d: EvalData<'a>) -> std::result::Result<$type_name, NameError> {
            if let $eval_data_name = d {
                Ok($v)
            } else {
                Err(NameError::Type($str, vec![d.into()]))
            }
        }
    };
}

impl_expect_type!(expect_bool, bool, EvalData::Bool(v), v, "bool");
impl_expect_type!(expect_int, isize, EvalData::Int(v), v, "int");
//impl_expect_type!(expect_float, f64, EvalData::Float(v), v, "float");
impl_expect_type!(expect_bytes, Cow<'a, [u8]>, EvalData::Bytes(v), v, "bytes");

fn map_eval_range<T>(
    range: &(Bound<Expr>, Bound<Expr>),
    read: &Read,
    use_qual: bool,
    mut map: impl FnMut(isize) -> T,
) -> std::result::Result<(Bound<T>, Bound<T>), NameError> {
    let start = match range.start_bound() {
        Bound::Included(s) => Bound::Included(map(expect_int(s.eval(read, use_qual)?)?)),
        Bound::Excluded(s) => Bound::Excluded(map(expect_int(s.eval(read, use_qual)?)?)),
        Bound::Unbounded => Bound::Unbounded,
    };
    let end = match range.end_bound() {
        Bound::Included(e) => Bound::Included(map(expect_int(e.eval(read, use_qual)?)?)),
        Bound::Excluded(e) => Bound::Excluded(map(expect_int(e.eval(read, use_qual)?)?)),
        Bound::Unbounded => Bound::Unbounded,
    };
    Ok((start, end))
}

fn range_inclusive_exclusive(range: &(Bound<usize>, Bound<usize>), len: usize) -> (usize, usize) {
    let start = match range.start_bound() {
        Bound::Included(&s) => s,
        Bound::Excluded(&s) => s + 1,
        Bound::Unbounded => 0,
    };
    let end = match range.end_bound() {
        Bound::Included(&e) => e + 1,
        Bound::Excluded(&e) => e,
        Bound::Unbounded => len,
    };
    (start, end)
}
