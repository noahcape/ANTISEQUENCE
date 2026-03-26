use rustc_hash::FxHashMap;
use smallvec::SmallVec;

use serde::Serialize;

use std::fmt;
use std::sync::Arc;

use crate::errors::{self, Name, NameError};
use crate::inline_string::*;

/// A FASTQ record as (name, sequence, quality) byte slices.
pub type FastqRecord<'a> = (&'a [u8], &'a [u8], &'a [u8]);

pub use End::*;

/// Left or right end.
#[derive(Copy, Clone, PartialEq, Debug)]
pub enum End {
    Left,
    Right,
}

/// Valid types of strings.
///
/// `Name` or `Seq` refer to the corresponding line in a fastq record.
/// Each Read contains multiple different strings of different types.
///
/// Uses 1-indexed conventions, like Name(1) and Seq(1), to follow fastq file naming conventions.
#[derive(Debug, Clone, Copy, PartialEq, Hash)]
pub enum StrType {
    Name(u8),
    Seq(u8),
}

/// A fastq read.
///
/// This is the core data structure that is manipulated by other operations in this library.
/// Both fastq records for paired-end reads are stored in the same `Read`.
#[derive(Debug, Clone)]
pub struct Read {
    str_mappings: Vec<(StrType, StrMappings)>,
}

/// A string and its correspondings mappings.
///
/// This is typically used to represent a name or sequence from a fastq record.
#[derive(Debug, Clone)]
pub struct StrMappings {
    mappings: SmallVec<[Mapping; 4]>,
    string: Vec<u8>,
    qual: Option<Vec<u8>>,

    // tracks where this string came from
    origin: Arc<Origin>,
    idx: usize,
}

impl StrMappings {
    #[inline(always)]
    pub fn new(string: Vec<u8>, origin: Arc<Origin>, idx: usize) -> Self {
        let mut mappings: SmallVec<[Mapping; 4]> = SmallVec::new();
        mappings.push(Mapping::new_default(string.len()));
        Self {
            mappings,
            string,
            qual: None,
            origin,
            idx,
        }
    }

    #[inline(always)]
    pub fn new_with_qual(string: Vec<u8>, qual: Vec<u8>, origin: Arc<Origin>, idx: usize) -> Self {
        let mut mappings: SmallVec<[Mapping; 4]> = SmallVec::new();
        mappings.push(Mapping::new_default(string.len()));
        Self {
            mappings,
            string,
            qual: Some(qual),
            origin,
            idx,
        }
    }

    #[inline(always)]
    pub fn data(&self, label: InlineString, attr: InlineString) -> Option<&Data> {
        self.mapping(label).and_then(|m| m.data(attr))
    }

    #[inline(always)]
    pub fn data_mut(&mut self, label: InlineString, attr: InlineString) -> Option<&mut Data> {
        self.mapping_mut(label).map(|m| m.data_mut(attr))
    }

    #[inline(always)]
    pub fn mapping(&self, label: InlineString) -> Option<&Mapping> {
        // iterate in reverse since recently added labels likely to be at the end
        self.mappings.iter().rev().find(|m| m.label == label)
    }

    #[inline(always)]
    pub fn mapping_mut(&mut self, label: InlineString) -> Option<&mut Mapping> {
        self.mappings.iter_mut().rev().find(|m| m.label == label)
    }

    #[inline(always)]
    pub fn add_mapping(&mut self, label: Option<InlineString>, start: usize, len: usize) {
        let Some(label) = label else {
            return;
        };
        if let Some(m) = self.mapping_mut(label) {
            m.start = start;
            m.len = len;
            if let Some(d) = m.data.as_mut() {
                d.clear();
            }
        } else {
            self.mappings.push(Mapping::new(label, start, len));
        }
    }

    #[inline(always)]
    pub fn string(&self) -> &[u8] {
        &self.string
    }

    #[inline(always)]
    pub fn string_mut(&mut self) -> &mut Vec<u8> {
        &mut self.string
    }

    #[inline(always)]
    pub fn qual(&self) -> Option<&[u8]> {
        self.qual.as_deref()
    }

    #[inline(always)]
    pub fn qual_mut(&mut self) -> Option<&mut Vec<u8>> {
        self.qual.as_mut()
    }

    #[inline(always)]
    pub fn substring(&self, mapping: &Mapping) -> &[u8] {
        &self.string[mapping.start..mapping.start + mapping.len]
    }

    #[inline(always)]
    pub fn substring_qual(&self, mapping: &Mapping) -> Option<&[u8]> {
        self.qual
            .as_ref()
            .map(|q| &q[mapping.start..mapping.start + mapping.len])
    }

    pub fn cut(
        &mut self,
        label: InlineString,
        new_label1: Option<InlineString>,
        new_label2: Option<InlineString>,
        cut_idx: isize,
    ) -> Result<(), NameError> {
        let (start, len) = {
            let mapping = self
                .mapping(label)
                .ok_or(NameError::NotInRead(Name::Label(label)))?;
            (mapping.start, mapping.len)
        };

        if cut_idx >= 0 {
            let cut = (cut_idx as usize).min(len);
            self.add_mapping(new_label1, start, cut);
            self.add_mapping(new_label2, start + cut, len - cut);
        } else {
            let cut = ((-cut_idx) as usize).min(len);
            self.add_mapping(new_label1, start, len - cut);
            self.add_mapping(new_label2, start + len - cut, cut);
        }

        Ok(())
    }

    pub fn intersect(
        &mut self,
        label1: InlineString,
        label2: InlineString,
        new_label: Option<InlineString>,
    ) -> Result<(), NameError> {
        let mapping1 = self
            .mapping(label1)
            .ok_or(NameError::NotInRead(Name::Label(label1)))?;
        let mapping2 = self
            .mapping(label2)
            .ok_or(NameError::NotInRead(Name::Label(label2)))?;

        if let Some((start, len)) = mapping1.intersection_interval(mapping2) {
            self.add_mapping(new_label, start, len);
        }

        Ok(())
    }

    pub fn union(
        &mut self,
        label1: InlineString,
        label2: InlineString,
        new_label: Option<InlineString>,
    ) -> Result<(), NameError> {
        let mapping1 = self
            .mapping(label1)
            .ok_or(NameError::NotInRead(Name::Label(label1)))?;
        let mapping2 = self
            .mapping(label2)
            .ok_or(NameError::NotInRead(Name::Label(label2)))?;

        let (start, len) = mapping1.union_interval(mapping2);
        self.add_mapping(new_label, start, len);

        Ok(())
    }

    pub fn set(
        &mut self,
        label: InlineString,
        new_str: &[u8],
        new_qual: Option<&[u8]>,
    ) -> Result<(), NameError> {
        let prev = self
            .mapping(label)
            .ok_or(NameError::NotInRead(Name::Label(label)))?
            .clone();

        self.mappings.iter_mut().for_each(|m| {
            if m.label.bytes().all(|b| b == b'*') {
                if new_str.len() >= prev.len {
                    m.len += new_str.len() - prev.len;
                } else {
                    m.len -= prev.len - new_str.len();
                }

                return;
            }

            use Intersection::*;
            match prev.intersect(m) {
                ABOverlap(len) => {
                    if len > new_str.len() {
                        m.start = prev.start;
                        m.len -= len - new_str.len();
                    } else if new_str.len() >= prev.len {
                        m.start += new_str.len() - prev.len;
                    } else {
                        m.start -= prev.len - new_str.len();
                    }
                }
                BAOverlap(len) => {
                    if len > new_str.len() {
                        m.len -= len - new_str.len();
                    }
                }
                AInsideB => {
                    if new_str.len() >= prev.len {
                        m.len += new_str.len() - prev.len;
                    } else {
                        m.len -= prev.len - new_str.len();
                    }
                }
                BInsideA => {
                    m.start = m.start.min(prev.start + new_str.len());
                    m.len = m.len.min(prev.start + new_str.len() - m.start);
                }
                Equal => {
                    m.len = new_str.len();
                }
                ABeforeB => {
                    if new_str.len() >= prev.len {
                        m.start += new_str.len() - prev.len;
                    } else {
                        m.start -= prev.len - new_str.len();
                    }
                }
                BBeforeA => (),
            }
        });

        // In-place replacement of bytes in self.string, minimizing reallocations.
        let start = prev.start;
        let old_len = prev.len;
        let new_len = new_str.len();
        // If new_str may alias self.string's buffer, copy into an owned temporary first to avoid UB
        let src_bytes: std::borrow::Cow<[u8]> = {
            let pr = self.string.as_ptr_range();
            let p = new_str.as_ptr();
            if p >= pr.start && p < pr.end {
                std::borrow::Cow::Owned(new_str.to_vec())
            } else {
                std::borrow::Cow::Borrowed(new_str)
            }
        };

        if new_len == old_len {
            if new_len != 0 {
                unsafe {
                    std::ptr::copy(
                        src_bytes.as_ref().as_ptr(),
                        self.string.as_mut_ptr().add(start),
                        new_len,
                    );
                }
            }
        } else if new_len < old_len {
            // write new bytes, then shift tail left and truncate
            if new_len != 0 {
                unsafe {
                    std::ptr::copy(
                        src_bytes.as_ref().as_ptr(),
                        self.string.as_mut_ptr().add(start),
                        new_len,
                    );
                }
            }
            let tail_src = start + old_len;
            let tail_dst = start + new_len;
            let tail_len = self.string.len() - tail_src;
            if tail_len > 0 {
                unsafe {
                    let base = self.string.as_mut_ptr();
                    std::ptr::copy(base.add(tail_src), base.add(tail_dst), tail_len);
                }
            }
            self.string
                .truncate(self.string.len() - (old_len - new_len));
        } else {
            // grow: make room by moving tail right, then write new bytes
            let diff = new_len - old_len;
            let tail_src = start + old_len;
            let tail_len = self.string.len() - tail_src;
            let orig_len = self.string.len();
            self.string.reserve(diff);
            self.string.resize(orig_len + diff, 0);
            if tail_len > 0 {
                unsafe {
                    let base = self.string.as_mut_ptr();
                    std::ptr::copy(base.add(tail_src), base.add(tail_src + diff), tail_len);
                }
            }
            if new_len != 0 {
                unsafe {
                    std::ptr::copy(
                        src_bytes.as_ref().as_ptr(),
                        self.string.as_mut_ptr().add(start),
                        new_len,
                    );
                }
            }
        }

        if let Some(qual_vec) = &mut self.qual {
            let qsrc = new_qual.expect("quality must be provided when qual is present");
            let q_new_len = qsrc.len();
            // alias check w.r.t qual buffer
            let q_bytes: std::borrow::Cow<[u8]> = {
                let pr = qual_vec.as_ptr_range();
                let p = qsrc.as_ptr();
                if p >= pr.start && p < pr.end {
                    std::borrow::Cow::Owned(qsrc.to_vec())
                } else {
                    std::borrow::Cow::Borrowed(qsrc)
                }
            };

            if q_new_len == old_len {
                if q_new_len != 0 {
                    unsafe {
                        std::ptr::copy(
                            q_bytes.as_ref().as_ptr(),
                            qual_vec.as_mut_ptr().add(start),
                            q_new_len,
                        );
                    }
                }
            } else if q_new_len < old_len {
                if q_new_len != 0 {
                    unsafe {
                        std::ptr::copy(
                            q_bytes.as_ref().as_ptr(),
                            qual_vec.as_mut_ptr().add(start),
                            q_new_len,
                        );
                    }
                }
                let tail_src = start + old_len;
                let tail_dst = start + q_new_len;
                let tail_len = qual_vec.len() - tail_src;
                if tail_len > 0 {
                    unsafe {
                        let base = qual_vec.as_mut_ptr();
                        std::ptr::copy(base.add(tail_src), base.add(tail_dst), tail_len);
                    }
                }
                qual_vec.truncate(qual_vec.len() - (old_len - q_new_len));
            } else {
                let diff = q_new_len - old_len;
                let tail_src = start + old_len;
                let tail_len = qual_vec.len() - tail_src;
                let orig_len = qual_vec.len();
                qual_vec.reserve(diff);
                qual_vec.resize(orig_len + diff, 0);
                if tail_len > 0 {
                    unsafe {
                        let base = qual_vec.as_mut_ptr();
                        std::ptr::copy(base.add(tail_src), base.add(tail_src + diff), tail_len);
                    }
                }
                if q_new_len != 0 {
                    unsafe {
                        std::ptr::copy(
                            q_bytes.as_ref().as_ptr(),
                            qual_vec.as_mut_ptr().add(start),
                            q_new_len,
                        );
                    }
                }
            }
        }

        Ok(())
    }

    pub fn trim(&mut self, label: InlineString) -> Result<(), NameError> {
        let trimmed = self
            .mapping(label)
            .ok_or(NameError::NotInRead(Name::Label(label)))?
            .clone();

        self.mappings.iter_mut().for_each(|m| {
            use Intersection::*;
            match trimmed.intersect(m) {
                ABOverlap(len) => {
                    m.start = trimmed.start;
                    m.len -= len;
                }
                BAOverlap(len) => {
                    m.len -= len;
                }
                AInsideB => {
                    m.len -= trimmed.len;
                }
                BInsideA => {
                    m.start = trimmed.start;
                    m.len = 0;
                }
                Equal => {
                    m.len = 0;
                }
                ABeforeB => {
                    m.start -= trimmed.len;
                }
                BBeforeA => (),
            }
        });

        // In-place remove [start, start+len) from string and qual with memmove + truncate
        let start = trimmed.start;
        let len = trimmed.len;
        let tail_src = start + len;
        let tail_len = self.string.len() - tail_src;
        if tail_len > 0 {
            unsafe {
                let base = self.string.as_mut_ptr();
                std::ptr::copy(base.add(tail_src), base.add(start), tail_len);
            }
        }
        self.string.truncate(self.string.len() - len);

        if let Some(qual_vec) = &mut self.qual {
            let tail_len_q = qual_vec.len() - tail_src;
            if tail_len_q > 0 {
                unsafe {
                    let base = qual_vec.as_mut_ptr();
                    std::ptr::copy(base.add(tail_src), base.add(start), tail_len_q);
                }
            }
            qual_vec.truncate(qual_vec.len() - len);
        }

        Ok(())
    }

    pub fn remove_internal(&mut self) {
        self.mappings
            .retain(|m| m.label.bytes().next() != Some(b'_'));
    }

    pub fn recycle(&mut self, string: Vec<u8>, origin: Arc<Origin>, idx: usize) {
        self.mappings.clear();
        self.mappings.push(Mapping::new_default(string.len()));
        self.string = string;
        self.qual = None;
        self.origin = origin;
        self.idx = idx;
    }

    pub fn recycle_with_qual(
        &mut self,
        string: Vec<u8>,
        qual: Vec<u8>,
        origin: Arc<Origin>,
        idx: usize,
    ) {
        self.mappings.clear();
        self.mappings.push(Mapping::new_default(string.len()));
        self.string = string;
        self.qual = Some(qual);
        self.origin = origin;
        self.idx = idx;
    }
}

/// A labeled mapping that corresponds to an interval/region in a string.
#[derive(Debug, Clone, PartialEq)]
pub struct Mapping {
    pub label: InlineString,
    pub start: usize,
    pub len: usize,
    data: Option<SmallAttrMap>,
}

/// Data types.
#[derive(Clone, PartialEq, Serialize)]
#[serde(untagged)]
pub enum Data {
    Bool(bool),
    Int(isize),
    Float(f64),
    Bytes(Vec<u8>),
}

#[derive(Debug, Clone, PartialEq)]
pub enum Intersection {
    ABOverlap(usize),
    BAOverlap(usize),
    AInsideB,
    BInsideA,
    Equal,
    ABeforeB,
    BBeforeA,
}

impl Mapping {
    #[inline(always)]
    pub fn new_default(len: usize) -> Self {
        Self {
            label: InlineString::new(b"*"),
            start: 0,
            len,
            data: None,
        }
    }

    #[inline(always)]
    pub fn new(label: InlineString, start: usize, len: usize) -> Self {
        Self {
            label,
            start,
            len,
            data: None,
        }
    }

    #[inline(always)]
    pub fn intersect(&self, b: &Self) -> Intersection {
        let a_start = self.start;
        let a_end = self.start + self.len;
        let b_start = b.start;
        let b_end = b.start + b.len;

        use Intersection::*;
        if a_start == b_start && a_end == b_end {
            Equal
        } else if a_start < b_start && b_end < a_end {
            BInsideA
        } else if b_start < a_start && a_end < b_end {
            AInsideB
        } else if a_start == b_start {
            if a_end > b_end {
                BAOverlap(b_end - a_start)
            } else {
                ABOverlap(a_end - b_start)
            }
        } else if a_end == b_end {
            if a_start > b_start {
                BAOverlap(b_end - a_start)
            } else {
                ABOverlap(a_end - b_start)
            }
        } else if a_start <= b_start && b_start < a_end {
            ABOverlap(a_end - b_start)
        } else if a_start < b_end && b_end <= a_end {
            BAOverlap(b_end - a_start)
        } else if a_end <= b_start {
            ABeforeB
        } else if b_end <= a_start {
            BBeforeA
        } else {
            unreachable!()
        }
    }

    #[inline(always)]
    pub fn intersection_interval(&self, b: &Self) -> Option<(usize, usize)> {
        let a_start = self.start;
        let a_end = self.start + self.len;
        let b_start = b.start;
        let b_end = b.start + b.len;

        if (b_start <= a_start && a_start < b_end) || (a_start <= b_start && b_start < a_end) {
            let start = a_start.max(b_start);
            let len = a_end.min(b_end) - start;
            Some((start, len))
        } else {
            None
        }
    }

    #[inline(always)]
    pub fn union_interval(&self, b: &Self) -> (usize, usize) {
        let a_start = self.start;
        let a_end = self.start + self.len;
        let b_start = b.start;
        let b_end = b.start + b.len;

        let start = a_start.min(b_start);
        let len = a_end.max(b_end) - start;
        (start, len)
    }

    #[inline(always)]
    pub fn data(&self, attr: InlineString) -> Option<&Data> {
        self.data.as_ref().and_then(|m| m.get(&attr))
    }

    #[inline(always)]
    pub fn data_mut(&mut self, attr: InlineString) -> &mut Data {
        self.data
            .get_or_insert_with(SmallAttrMap::default)
            .get_or_insert_default(attr)
    }

    pub fn remove_data(&mut self, attr: &InlineString) {
        if let Some(d) = &mut self.data {
            d.remove(attr);
        }
    }
}

#[derive(Debug, Clone, PartialEq)]
struct SmallAttrMap {
    small: SmallVec<[(InlineString, Data); 4]>,
    map: Option<FxHashMap<InlineString, Data>>,
}

impl Default for SmallAttrMap {
    fn default() -> Self {
        Self {
            small: SmallVec::new(),
            map: None,
        }
    }
}

impl SmallAttrMap {
    fn clear(&mut self) {
        self.small.clear();
        if let Some(m) = &mut self.map {
            m.clear();
        }
    }

    fn get(&self, attr: &InlineString) -> Option<&Data> {
        if let Some(m) = &self.map {
            return m.get(attr);
        }
        self.small
            .iter()
            .find_map(|(k, v)| if k == attr { Some(v) } else { None })
    }

    fn get_or_insert_default(&mut self, attr: InlineString) -> &mut Data {
        // If we already promoted to map, insert there via helper to avoid borrow conflicts.
        if self.map.is_some() {
            return self.entry_in_map(attr);
        }

        // Try to find in small via index to avoid overlapping borrows
        let len_now = self.small.len();
        for i in 0..len_now {
            if self.small[i].0 == attr {
                return &mut self.small[i].1;
            }
        }

        // If room left inline, push
        let inline_cap = self.small.inline_size();
        if self.small.len() < inline_cap {
            self.small.push((attr, Data::Bool(false)));
            let idx = self.small.len() - 1;
            return &mut self.small[idx].1;
        }

        // Promote to hashmap
        let mut m: FxHashMap<InlineString, Data> = FxHashMap::default();
        m.reserve(self.small.len() + 1);
        for (k, v) in self.small.drain(..) {
            m.insert(k, v);
        }
        self.map = Some(m);
        // Now safely insert and return from the map
        self.entry_in_map(attr)
    }

    fn remove(&mut self, attr: &InlineString) {
        if let Some(m) = &mut self.map {
            m.remove(attr);
        } else {
            self.small.retain(|(k, _)| k != attr);
        }
    }

    fn for_each<F: FnMut(&InlineString, &Data)>(&self, mut f: F) {
        if let Some(m) = &self.map {
            for (k, v) in m.iter() {
                f(k, v);
            }
        } else {
            for (k, v) in self.small.iter() {
                f(k, v);
            }
        }
    }

    #[inline(always)]
    fn entry_in_map(&mut self, attr: InlineString) -> &mut Data {
        self.map
            .as_mut()
            .unwrap()
            .entry(attr)
            .or_insert(Data::Bool(false))
    }
}

impl Default for Read {
    fn default() -> Self {
        Self::new()
    }
}

impl Read {
    #[inline(always)]
    pub fn new() -> Self {
        Self {
            str_mappings: Vec::with_capacity(4),
        }
    }

    #[inline(always)]
    pub fn clear(&mut self) {
        self.str_mappings.clear();
    }

    #[inline(always)]
    pub fn has_names(&self, names: &[crate::expr::LabelOrAttr]) -> bool {
        for name in names {
            match name {
                crate::expr::LabelOrAttr::Label(l) => {
                    if let Some(s) = self.str_mappings(l.str_type) {
                        if s.mapping(l.label).is_some() {
                            continue;
                        }
                    }
                    return false;
                }
                crate::expr::LabelOrAttr::Attr(a) => {
                    if let Some(s) = self.str_mappings(a.str_type) {
                        if s.data(a.label, a.attr).is_some() {
                            continue;
                        }
                    }
                    return false;
                }
            }
        }

        true
    }

    pub fn add_fastq(
        &mut self,
        str_type_idx: u8,
        name: &[u8],
        seq: &[u8],
        qual: &[u8],
        origin: Arc<Origin>,
        idx: usize,
    ) {
        let name = StrMappings::new(name.to_owned(), Arc::clone(&origin), idx);
        let seq = StrMappings::new_with_qual(seq.to_owned(), qual.to_owned(), origin, idx);
        self.str_mappings.push((StrType::Name(str_type_idx), name));
        self.str_mappings.push((StrType::Seq(str_type_idx), seq));
    }

    pub fn add_fastq_parts(
        &mut self,
        str_type_idx: u8,
        name: Option<&[u8]>,
        seq: &[u8],
        qual: Option<&[u8]>,
        origin: Arc<Origin>,
        idx: usize,
    ) {
        if let Some(n) = name {
            let name_sm = StrMappings::new(n.to_owned(), Arc::clone(&origin), idx);
            self.str_mappings
                .push((StrType::Name(str_type_idx), name_sm));
        }
        let seq_sm = match qual {
            Some(q) => StrMappings::new_with_qual(seq.to_owned(), q.to_owned(), origin, idx),
            None => StrMappings::new(seq.to_owned(), origin, idx),
        };
        self.str_mappings.push((StrType::Seq(str_type_idx), seq_sm));
    }

    pub fn add_fastq_recycled(
        &mut self,
        str_type_idx: u8,
        name: &[u8],
        seq: &[u8],
        qual: &[u8],
        origin: Arc<Origin>,
        idx: usize,
    ) {
        // Try to reuse existing StrMappings if they exist
        let name_type = StrType::Name(str_type_idx);
        let seq_type = StrType::Seq(str_type_idx);

        let mut name_found = false;
        let mut seq_found = false;

        for (t, sm) in &mut self.str_mappings {
            if *t == name_type {
                let mut s = std::mem::take(&mut sm.string);
                s.clear();
                s.extend_from_slice(name);
                sm.recycle(s, Arc::clone(&origin), idx);
                name_found = true;
            } else if *t == seq_type {
                let mut s = std::mem::take(&mut sm.string);
                s.clear();
                s.extend_from_slice(seq);

                let mut q = sm.qual.take().unwrap_or_default();
                q.clear();
                q.extend_from_slice(qual);

                sm.recycle_with_qual(s, q, Arc::clone(&origin), idx);
                seq_found = true;
            }
        }

        if !name_found {
            let name = StrMappings::new(name.to_owned(), Arc::clone(&origin), idx);
            self.str_mappings.push((name_type, name));
        }
        if !seq_found {
            let seq = StrMappings::new_with_qual(seq.to_owned(), qual.to_owned(), origin, idx);
            self.str_mappings.push((seq_type, seq));
        }
    }

    pub fn add_fastq_parts_recycled(
        &mut self,
        str_type_idx: u8,
        name: Option<&[u8]>,
        seq: &[u8],
        qual: Option<&[u8]>,
        origin: Arc<Origin>,
        idx: usize,
    ) {
        let name_type = StrType::Name(str_type_idx);
        let seq_type = StrType::Seq(str_type_idx);

        let mut name_found = false;
        let mut seq_found = false;

        for (t, sm) in &mut self.str_mappings {
            if *t == name_type {
                if let Some(n) = name {
                    let mut s = std::mem::take(&mut sm.string);
                    s.clear();
                    s.extend_from_slice(n);
                    sm.recycle(s, Arc::clone(&origin), idx);
                    name_found = true;
                }
            } else if *t == seq_type {
                let mut s = std::mem::take(&mut sm.string);
                s.clear();
                s.extend_from_slice(seq);

                if let Some(q_bytes) = qual {
                    let mut q = sm.qual.take().unwrap_or_default();
                    q.clear();
                    q.extend_from_slice(q_bytes);
                    sm.recycle_with_qual(s, q, Arc::clone(&origin), idx);
                } else {
                    sm.recycle(s, Arc::clone(&origin), idx);
                }
                seq_found = true;
            }
        }

        if !name_found {
            if let Some(n) = name {
                let name_sm = StrMappings::new(n.to_owned(), Arc::clone(&origin), idx);
                self.str_mappings.push((name_type, name_sm));
            }
        }
        if !seq_found {
            let seq_sm = match qual {
                Some(q) => StrMappings::new_with_qual(seq.to_owned(), q.to_owned(), origin, idx),
                None => StrMappings::new(seq.to_owned(), origin, idx),
            };
            self.str_mappings.push((seq_type, seq_sm));
        }
    }

    pub fn set_fastq_entry(
        &mut self,
        slot_idx: usize,
        str_type: StrType,
        string: &[u8],
        qual: Option<&[u8]>,
        origin: Arc<Origin>,
        idx: usize,
    ) {
        if slot_idx < self.str_mappings.len() {
            let (t, sm) = &mut self.str_mappings[slot_idx];
            *t = str_type;

            let mut s = std::mem::take(&mut sm.string);
            s.clear();
            s.extend_from_slice(string);

            if let Some(q_bytes) = qual {
                let mut q = sm.qual.take().unwrap_or_default();
                q.clear();
                q.extend_from_slice(q_bytes);
                sm.recycle_with_qual(s, q, origin, idx);
            } else {
                sm.recycle(s, origin, idx);
            }
        } else {
            let sm = if let Some(q) = qual {
                StrMappings::new_with_qual(string.to_owned(), q.to_owned(), origin, idx)
            } else {
                StrMappings::new(string.to_owned(), origin, idx)
            };
            self.str_mappings.push((str_type, sm));
        }
    }

    /// Returns (name, sequence, quality) for the given string type index.
    #[inline(always)]
    pub fn to_fastq(&self, str_type_idx: u8) -> Result<FastqRecord<'_>, NameError> {
        let name = self
            .str_mappings(StrType::Name(str_type_idx))
            .ok_or(NameError::NotInRead(Name::StrType(StrType::Name(
                str_type_idx,
            ))))?;
        let seq = self.str_mappings(StrType::Seq(str_type_idx)).unwrap();
        Ok((name.string(), seq.string(), seq.qual().unwrap()))
    }

    pub fn to_json(&self) -> String {
        serde_json::to_string(&SerializableRead::from(self)).unwrap()
    }

    #[inline(always)]
    pub fn str_mappings(&self, str_type: StrType) -> Option<&StrMappings> {
        self.str_mappings
            .iter()
            .find_map(|(t, m)| if *t == str_type { Some(m) } else { None })
    }

    #[inline(always)]
    pub fn str_mappings_mut(&mut self, str_type: StrType) -> Option<&mut StrMappings> {
        self.str_mappings
            .iter_mut()
            .find_map(|(t, m)| if *t == str_type { Some(m) } else { None })
    }

    #[inline(always)]
    pub fn mapping(&self, str_type: StrType, label: InlineString) -> Result<&Mapping, NameError> {
        self.str_mappings(str_type)
            .ok_or(NameError::NotInRead(Name::StrType(str_type)))?
            .mapping(label)
            .ok_or(NameError::NotInRead(Name::Label(label)))
    }

    #[inline(always)]
    pub fn mapping_mut(
        &mut self,
        str_type: StrType,
        label: InlineString,
    ) -> Result<&mut Mapping, NameError> {
        self.str_mappings_mut(str_type)
            .ok_or(NameError::NotInRead(Name::StrType(str_type)))?
            .mapping_mut(label)
            .ok_or(NameError::NotInRead(Name::Label(label)))
    }

    pub fn data(
        &self,
        str_type: StrType,
        label: InlineString,
        attr: InlineString,
    ) -> Result<&Data, NameError> {
        self.str_mappings(str_type)
            .ok_or(NameError::NotInRead(Name::StrType(str_type)))?
            .mapping(label)
            .ok_or(NameError::NotInRead(Name::Label(label)))?
            .data(attr)
            .ok_or(NameError::NotInRead(Name::Attr(attr)))
    }

    pub fn data_mut(
        &mut self,
        str_type: StrType,
        label: InlineString,
        attr: InlineString,
    ) -> Result<&mut Data, NameError> {
        Ok(self
            .str_mappings_mut(str_type)
            .ok_or(NameError::NotInRead(Name::StrType(str_type)))?
            .mapping_mut(label)
            .ok_or(NameError::NotInRead(Name::Label(label)))?
            .data_mut(attr))
    }

    pub fn remove_data(&mut self, str_type: StrType, label: InlineString, attr: &InlineString) {
        if let Some(sm) = self.str_mappings_mut(str_type) {
            if let Some(m) = sm.mapping_mut(label) {
                m.remove_data(attr);
            }
        }
    }

    #[inline(always)]
    pub fn substring(&self, str_type: StrType, label: InlineString) -> Result<&[u8], NameError> {
        let str_mappings = self
            .str_mappings(str_type)
            .ok_or(NameError::NotInRead(Name::StrType(str_type)))?;
        let mapping = str_mappings
            .mapping(label)
            .ok_or(NameError::NotInRead(Name::Label(label)))?;
        Ok(str_mappings.substring(mapping))
    }

    #[inline(always)]
    pub fn substring_qual(
        &self,
        str_type: StrType,
        label: InlineString,
    ) -> Result<Option<&[u8]>, NameError> {
        let str_mappings = self
            .str_mappings(str_type)
            .ok_or(NameError::NotInRead(Name::StrType(str_type)))?;
        let mapping = str_mappings
            .mapping(label)
            .ok_or(NameError::NotInRead(Name::Label(label)))?;
        Ok(str_mappings.substring_qual(mapping))
    }

    pub fn cut(
        &mut self,
        str_type: StrType,
        label: InlineString,
        new_label1: Option<InlineString>,
        new_label2: Option<InlineString>,
        cut_idx: isize,
    ) -> Result<(), NameError> {
        self.str_mappings_mut(str_type)
            .ok_or(NameError::NotInRead(Name::StrType(str_type)))?
            .cut(label, new_label1, new_label2, cut_idx)
    }

    pub fn intersect(
        &mut self,
        str_type: StrType,
        label1: InlineString,
        label2: InlineString,
        new_label: Option<InlineString>,
    ) -> Result<(), NameError> {
        self.str_mappings_mut(str_type)
            .ok_or(NameError::NotInRead(Name::StrType(str_type)))?
            .intersect(label1, label2, new_label)
    }

    pub fn union(
        &mut self,
        str_type: StrType,
        label1: InlineString,
        label2: InlineString,
        new_label: Option<InlineString>,
    ) -> Result<(), NameError> {
        self.str_mappings_mut(str_type)
            .ok_or(NameError::NotInRead(Name::StrType(str_type)))?
            .union(label1, label2, new_label)
    }

    pub fn set(
        &mut self,
        str_type: StrType,
        label: InlineString,
        new_str: &[u8],
        new_qual: Option<&[u8]>,
    ) -> Result<(), NameError> {
        self.str_mappings_mut(str_type)
            .ok_or(NameError::NotInRead(Name::StrType(str_type)))?
            .set(label, new_str, new_qual)
    }

    pub fn trim(&mut self, str_type: StrType, label: InlineString) -> Result<(), NameError> {
        self.str_mappings_mut(str_type)
            .ok_or(NameError::NotInRead(Name::StrType(str_type)))?
            .trim(label)
    }

    pub fn remove_internal(&mut self) {
        self.str_mappings
            .iter_mut()
            .for_each(|(_, s)| s.remove_internal());
    }

    pub fn first_idx(&self) -> usize {
        self.str_mappings.iter().map(|(_, s)| s.idx).min().unwrap()
    }
}

impl Data {
    pub fn as_bool(&self) -> bool {
        use Data::*;
        match self {
            Bool(x) => *x,
            Int(x) => *x != 0,
            Float(x) => *x != 0.0,
            Bytes(x) => !x.is_empty(),
        }
    }

    pub fn as_int(&self) -> Result<isize, NameError> {
        use Data::*;
        match self {
            Bool(x) => Ok(if *x { 1 } else { 0 }),
            Int(x) => Ok(*x),
            Float(x) => Ok(*x as isize),
            Bytes(_) => Err(NameError::Type("bool or uint", vec![self.clone()])),
        }
    }

    pub fn len(&self) -> Result<usize, NameError> {
        use Data::*;
        match self {
            Bool(_) => Err(NameError::Type("bytes", vec![self.clone()])),
            Int(_) => Err(NameError::Type("bytes", vec![self.clone()])),
            Float(_) => Err(NameError::Type("bytes", vec![self.clone()])),
            Bytes(x) => Ok(x.len()),
        }
    }

    pub fn is_empty(&self) -> Result<bool, NameError> {
        self.len().map(|l| l == 0)
    }
}

impl fmt::Display for StrMappings {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        use colored::Colorize;

        let len = self
            .mappings
            .iter()
            .map(|m| m.label.len())
            .max()
            .unwrap()
            .max(5);

        for m in &self.mappings {
            let curr = if m.len == 0 {
                let mut c = vec![b' '; self.string().len() + 1];
                c[m.start] = b'.';
                String::from_utf8(c).unwrap()
            } else {
                let mut c = vec![b' '; self.string().len() + 1];
                c[m.start..m.start + m.len].fill(b'-');
                c[m.start] = b'|';
                c[m.start + m.len - 1] = b'|';
                String::from_utf8(c).unwrap()
            };

            if m.label.bytes().next() == Some(b'_') {
                // internal
                write!(
                    f,
                    " {: <len$} {}",
                    m.label.to_string().bold().dimmed(),
                    curr.dimmed()
                )?;
            } else {
                write!(f, " {: <len$} {}", m.label.to_string().bold(), curr)?;
            }

            if let Some(data) = &m.data {
                data.for_each(|k, v| {
                    let _ = write!(f, " {}={}", k.to_string().bold(), v);
                });
            }
            writeln!(f)?;
        }

        writeln!(
            f,
            " {: <len$} {}",
            "str:".bold().green(),
            std::str::from_utf8(self.string()).unwrap().green()
        )?;

        if let Some(qual) = self.qual() {
            writeln!(
                f,
                " {: <len$} {}",
                "qual:".bold().green(),
                std::str::from_utf8(qual).unwrap().green()
            )?;
        }

        write!(
            f,
            " {: <len$} record {} in {}",
            "from:".bold(),
            self.idx,
            &*self.origin
        )?;

        Ok(())
    }
}

impl fmt::Display for Read {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        use colored::Colorize;

        for (str_type, str_mapping) in &self.str_mappings {
            writeln!(
                f,
                "{}:\n{}",
                str_type.to_string().bold().underline(),
                str_mapping
            )?;
        }
        Ok(())
    }
}

impl fmt::Display for Data {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        use Data::*;
        match self {
            Bool(x) => write!(f, "{}", x),
            Int(x) => write!(f, "{}", x),
            Float(x) => write!(f, "{}", x),
            Bytes(x) => write!(f, "{}", std::str::from_utf8(x).unwrap()),
        }
    }
}

impl fmt::Debug for Data {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        use Data::*;
        match self {
            Bool(x) => write!(f, "bool {}", x),
            Int(x) => write!(f, "int {}", x),
            Float(x) => write!(f, "float {}", x),
            Bytes(x) => write!(f, "bytes \"{}\"", std::str::from_utf8(x).unwrap()),
        }
    }
}

impl StrType {
    pub fn new(str_type: &[u8]) -> Result<Self, errors::Error> {
        use StrType::*;
        if str_type.starts_with(b"name") {
            let i = std::str::from_utf8(&str_type[4..])
                .unwrap()
                .parse::<u8>()
                .map_err(|_| errors::Error::Parse {
                    string: errors::utf8(str_type),
                    context: errors::utf8(str_type),
                    reason: "not a known valid string type. Expected \"name1\", \"name2\", etc.",
                })?;
            Ok(Name(i))
        } else if str_type.starts_with(b"seq") {
            let i = std::str::from_utf8(&str_type[3..])
                .unwrap()
                .parse::<u8>()
                .map_err(|_| errors::Error::Parse {
                    string: errors::utf8(str_type),
                    context: errors::utf8(str_type),
                    reason: "not a known valid string type. Expected \"seq1\", \"seq2\", etc.",
                })?;
            Ok(Seq(i))
        } else {
            Err(errors::Error::Parse {
                string: errors::utf8(str_type),
                context: errors::utf8(str_type),
                reason: "not a known valid string type. Expected \"name1\", \"seq1\", etc.",
            })
        }
    }
}

impl fmt::Display for StrType {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        use StrType::*;
        match self {
            Name(i) => write!(f, "name{i}"),
            Seq(i) => write!(f, "seq{i}"),
        }
    }
}

#[derive(Debug, Clone)]
pub enum Origin {
    File(String),
    Bytes,
}

impl fmt::Display for Origin {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match self {
            Origin::File(file) => write!(f, "file: \"{}\"", file),
            Origin::Bytes => write!(f, "bytes"),
        }
    }
}

impl End {
    pub fn to_isize(&self, i: usize) -> isize {
        match self {
            Left => i as isize,
            Right => -(i as isize),
        }
    }
}

#[derive(Serialize)]
pub(crate) struct SerializableRead(FxHashMap<String, SerializableStrMapping>);

#[derive(Serialize)]
struct SerializableStrMapping {
    #[serde(flatten)]
    mappings: FxHashMap<String, SerializableMapping>,
    idx: usize,
}

#[derive(Serialize)]
struct SerializableMapping {
    string: String,
    #[serde(skip_serializing_if = "Option::is_none")]
    qual: Option<String>,
    #[serde(skip_serializing_if = "FxHashMap::is_empty")]
    data: FxHashMap<String, Data>,
}

impl From<&Read> for SerializableRead {
    fn from(r: &Read) -> Self {
        let mut str_mappings = FxHashMap::default();

        for (str_type, str_mapping) in &r.str_mappings {
            let mut mappings = FxHashMap::default();

            for mapping in &str_mapping.mappings {
                if mapping.label.bytes().next() == Some(b'_') {
                    continue;
                }

                let data = mapping
                    .data
                    .as_ref()
                    .map(|m| {
                        let mut out = FxHashMap::default();
                        m.for_each(|attr, value| {
                            out.insert(attr.to_string(), value.clone());
                        });
                        out
                    })
                    .unwrap_or_default();
                let serializable_mapping = SerializableMapping {
                    string: std::str::from_utf8(str_mapping.substring(mapping))
                        .unwrap()
                        .to_string(),
                    qual: str_mapping
                        .substring_qual(mapping)
                        .map(|q| std::str::from_utf8(q).unwrap().to_string()),
                    data,
                };
                mappings.insert(mapping.label.to_string(), serializable_mapping);
            }

            let serializable_str_mapping = SerializableStrMapping {
                mappings,
                idx: str_mapping.idx,
            };
            str_mappings.insert(str_type.to_string(), serializable_str_mapping);
        }

        Self(str_mappings)
    }
}

#[cfg(test)]
mod size_tests {
    use super::Read;

    #[test]
    fn read_size_guard() {
        let sz = std::mem::size_of::<Read>();
        assert!(sz <= 128, "Read is too large: {} bytes", sz);
    }
}

#[cfg(test)]
mod read_tests {
    use super::*;
    use std::sync::Arc;

    fn test_origin() -> Arc<Origin> {
        Arc::new(Origin::File("test.fastq".to_string()))
    }

    #[test]
    fn test_str_mappings_new() {
        let origin = test_origin();
        let sm = StrMappings::new(b"ACGT".to_vec(), origin, 0);
        assert_eq!(sm.string(), b"ACGT");
        assert!(sm.qual().is_none());
        assert!(sm.mapping(InlineString::new(b"*")).is_some());
    }

    #[test]
    fn test_str_mappings_new_with_qual() {
        let origin = test_origin();
        let sm = StrMappings::new_with_qual(b"ACGT".to_vec(), b"IIII".to_vec(), origin, 0);
        assert_eq!(sm.string(), b"ACGT");
        assert_eq!(sm.qual(), Some(b"IIII".as_slice()));
    }

    #[test]
    fn test_str_mappings_add_mapping() {
        let origin = test_origin();
        let mut sm = StrMappings::new(b"ACGTACGT".to_vec(), origin, 0);
        sm.add_mapping(Some(InlineString::new(b"test")), 2, 4);
        let mapping = sm.mapping(InlineString::new(b"test")).unwrap();
        assert_eq!(mapping.start, 2);
        assert_eq!(mapping.len, 4);
        assert_eq!(sm.substring(mapping), b"GTAC");
    }

    #[test]
    fn test_str_mappings_add_mapping_none() {
        let origin = test_origin();
        let mut sm = StrMappings::new(b"ACGT".to_vec(), origin, 0);
        let count_before = sm.mappings.len();
        sm.add_mapping(None, 0, 2);
        assert_eq!(sm.mappings.len(), count_before);
    }

    #[test]
    fn test_str_mappings_substring_qual() {
        let origin = test_origin();
        let mut sm =
            StrMappings::new_with_qual(b"ACGTACGT".to_vec(), b"IIIIIIII".to_vec(), origin, 0);
        sm.add_mapping(Some(InlineString::new(b"test")), 2, 4);
        let mapping = sm.mapping(InlineString::new(b"test")).unwrap();
        assert_eq!(sm.substring_qual(mapping), Some(b"IIII".as_slice()));
    }

    #[test]
    fn test_str_mappings_cut_positive() {
        let origin = test_origin();
        let mut sm = StrMappings::new(b"ACGTACGT".to_vec(), origin, 0);
        sm.cut(
            InlineString::new(b"*"),
            Some(InlineString::new(b"left")),
            Some(InlineString::new(b"right")),
            3,
        )
        .unwrap();
        let left = sm.mapping(InlineString::new(b"left")).unwrap();
        let right = sm.mapping(InlineString::new(b"right")).unwrap();
        assert_eq!(sm.substring(left), b"ACG");
        assert_eq!(sm.substring(right), b"TACGT");
    }

    #[test]
    fn test_str_mappings_cut_negative() {
        let origin = test_origin();
        let mut sm = StrMappings::new(b"ACGTACGT".to_vec(), origin, 0);
        sm.cut(
            InlineString::new(b"*"),
            Some(InlineString::new(b"left")),
            Some(InlineString::new(b"right")),
            -3,
        )
        .unwrap();
        let left = sm.mapping(InlineString::new(b"left")).unwrap();
        let right = sm.mapping(InlineString::new(b"right")).unwrap();
        assert_eq!(sm.substring(left), b"ACGTA");
        assert_eq!(sm.substring(right), b"CGT");
    }

    #[test]
    fn test_str_mappings_cut_missing_label() {
        let origin = test_origin();
        let mut sm = StrMappings::new(b"ACGT".to_vec(), origin, 0);
        let result = sm.cut(
            InlineString::new(b"nonexistent"),
            Some(InlineString::new(b"left")),
            Some(InlineString::new(b"right")),
            2,
        );
        assert!(result.is_err());
    }

    #[test]
    fn test_str_mappings_intersect() {
        let origin = test_origin();
        let mut sm = StrMappings::new(b"ACGTACGT".to_vec(), origin, 0);
        sm.add_mapping(Some(InlineString::new(b"a")), 0, 5);
        sm.add_mapping(Some(InlineString::new(b"b")), 3, 5);
        sm.intersect(
            InlineString::new(b"a"),
            InlineString::new(b"b"),
            Some(InlineString::new(b"inter")),
        )
        .unwrap();
        let inter = sm.mapping(InlineString::new(b"inter")).unwrap();
        assert_eq!(inter.start, 3);
        assert_eq!(inter.len, 2);
    }

    #[test]
    fn test_str_mappings_union() {
        let origin = test_origin();
        let mut sm = StrMappings::new(b"ACGTACGT".to_vec(), origin, 0);
        sm.add_mapping(Some(InlineString::new(b"a")), 0, 3);
        sm.add_mapping(Some(InlineString::new(b"b")), 5, 3);
        sm.union(
            InlineString::new(b"a"),
            InlineString::new(b"b"),
            Some(InlineString::new(b"uni")),
        )
        .unwrap();
        let uni = sm.mapping(InlineString::new(b"uni")).unwrap();
        assert_eq!(uni.start, 0);
        assert_eq!(uni.len, 8);
    }

    #[test]
    fn test_mapping_new() {
        let m = Mapping::new(InlineString::new(b"test"), 5, 10);
        assert_eq!(m.start, 5);
        assert_eq!(m.len, 10);
    }

    #[test]
    fn test_mapping_new_default() {
        let m = Mapping::new_default(20);
        assert_eq!(m.start, 0);
        assert_eq!(m.len, 20);
    }

    #[test]
    fn test_mapping_intersection_interval() {
        let m1 = Mapping::new(InlineString::new(b"a"), 0, 5);
        let m2 = Mapping::new(InlineString::new(b"b"), 3, 5);
        let (start, len) = m1.intersection_interval(&m2).unwrap();
        assert_eq!(start, 3);
        assert_eq!(len, 2);
    }

    #[test]
    fn test_mapping_intersection_interval_no_overlap() {
        let m1 = Mapping::new(InlineString::new(b"a"), 0, 3);
        let m2 = Mapping::new(InlineString::new(b"b"), 5, 3);
        assert!(m1.intersection_interval(&m2).is_none());
    }

    #[test]
    fn test_mapping_union_interval() {
        let m1 = Mapping::new(InlineString::new(b"a"), 0, 3);
        let m2 = Mapping::new(InlineString::new(b"b"), 5, 3);
        let (start, len) = m1.union_interval(&m2);
        assert_eq!(start, 0);
        assert_eq!(len, 8);
    }

    #[test]
    fn test_read_new() {
        let read = Read::new();
        assert!(read.str_mappings.is_empty());
    }

    #[test]
    fn test_read_add_fastq() {
        let origin = test_origin();
        let mut read = Read::new();
        read.add_fastq(1, b"read1", b"ACGT", b"IIII", origin, 0);
        assert!(read.str_mappings(StrType::Seq(1)).is_some());
        assert!(read.str_mappings(StrType::Name(1)).is_some());
    }

    #[test]
    fn test_str_type_display() {
        assert_eq!(format!("{}", StrType::Seq(1)), "seq1");
        assert_eq!(format!("{}", StrType::Seq(2)), "seq2");
        assert_eq!(format!("{}", StrType::Name(1)), "name1");
    }

    #[test]
    fn test_str_type_new() {
        assert!(StrType::new(b"seq1").is_ok());
        assert!(StrType::new(b"seq2").is_ok());
        assert!(StrType::new(b"name1").is_ok());
        assert!(StrType::new(b"invalid").is_err());
    }

    #[test]
    fn test_end_enum() {
        assert_eq!(Left, End::Left);
        assert_eq!(Right, End::Right);
    }

    #[test]
    fn test_end_to_isize() {
        assert_eq!(Left.to_isize(5), 5);
        assert_eq!(Right.to_isize(5), -5);
    }

    #[test]
    fn test_origin_display() {
        let origin = Origin::File("test.fastq".to_string());
        assert!(format!("{}", origin).contains("test.fastq"));
        let origin_bytes = Origin::Bytes;
        assert_eq!(format!("{}", origin_bytes), "bytes");
    }

    #[test]
    fn test_intersection_ab_overlap() {
        let m1 = Mapping::new(InlineString::new(b"a"), 0, 5);
        let m2 = Mapping::new(InlineString::new(b"b"), 3, 5);
        let inter = m1.intersect(&m2);
        assert!(matches!(inter, Intersection::ABOverlap(_)));
    }

    #[test]
    fn test_intersection_ba_overlap() {
        let m1 = Mapping::new(InlineString::new(b"a"), 3, 5);
        let m2 = Mapping::new(InlineString::new(b"b"), 0, 5);
        let inter = m1.intersect(&m2);
        assert!(matches!(inter, Intersection::BAOverlap(_)));
    }

    #[test]
    fn test_intersection_a_before_b() {
        let m1 = Mapping::new(InlineString::new(b"a"), 0, 3);
        let m2 = Mapping::new(InlineString::new(b"b"), 5, 3);
        let inter = m1.intersect(&m2);
        assert!(matches!(inter, Intersection::ABeforeB));
    }

    #[test]
    fn test_intersection_b_before_a() {
        let m1 = Mapping::new(InlineString::new(b"a"), 5, 3);
        let m2 = Mapping::new(InlineString::new(b"b"), 0, 3);
        let inter = m1.intersect(&m2);
        assert!(matches!(inter, Intersection::BBeforeA));
    }

    #[test]
    fn test_intersection_a_inside_b() {
        let m1 = Mapping::new(InlineString::new(b"a"), 2, 3);
        let m2 = Mapping::new(InlineString::new(b"b"), 0, 10);
        let inter = m1.intersect(&m2);
        assert!(matches!(inter, Intersection::AInsideB));
    }

    #[test]
    fn test_intersection_b_inside_a() {
        let m1 = Mapping::new(InlineString::new(b"a"), 0, 10);
        let m2 = Mapping::new(InlineString::new(b"b"), 2, 3);
        let inter = m1.intersect(&m2);
        assert!(matches!(inter, Intersection::BInsideA));
    }

    #[test]
    fn test_intersection_equal() {
        let m1 = Mapping::new(InlineString::new(b"a"), 0, 5);
        let m2 = Mapping::new(InlineString::new(b"b"), 0, 5);
        let inter = m1.intersect(&m2);
        assert!(matches!(inter, Intersection::Equal));
    }

    // -- Mapping data methods --

    #[test]
    fn test_mapping_data_none() {
        let m = Mapping::new(InlineString::new(b"test"), 0, 5);
        assert!(m.data(InlineString::new(b"attr")).is_none());
    }

    #[test]
    fn test_mapping_data_mut() {
        let mut m = Mapping::new(InlineString::new(b"test"), 0, 5);
        *m.data_mut(InlineString::new(b"score")) = Data::Int(42);
        assert_eq!(m.data(InlineString::new(b"score")), Some(&Data::Int(42)));
    }

    #[test]
    fn test_mapping_data_mut_multiple() {
        let mut m = Mapping::new(InlineString::new(b"test"), 0, 5);
        *m.data_mut(InlineString::new(b"a")) = Data::Int(1);
        *m.data_mut(InlineString::new(b"b")) = Data::Int(2);
        *m.data_mut(InlineString::new(b"c")) = Data::Int(3);
        *m.data_mut(InlineString::new(b"d")) = Data::Int(4);
        // This should trigger promotion to hashmap (inline_size = 4)
        *m.data_mut(InlineString::new(b"e")) = Data::Int(5);
        assert_eq!(m.data(InlineString::new(b"a")), Some(&Data::Int(1)));
        assert_eq!(m.data(InlineString::new(b"e")), Some(&Data::Int(5)));
    }

    // -- Data methods --

    #[test]
    fn test_data_as_bool() {
        assert!(Data::Bool(true).as_bool());
        assert!(!(Data::Bool(false).as_bool()));
        assert!(Data::Int(1).as_bool());
        assert!(!(Data::Int(0).as_bool()));
        assert!(Data::Float(1.0).as_bool());
        assert!(!(Data::Float(0.0).as_bool()));
        assert!(Data::Bytes(b"test".to_vec()).as_bool());
        assert!(!(Data::Bytes(vec![]).as_bool()));
    }

    #[test]
    fn test_data_as_int() {
        assert_eq!(Data::Bool(true).as_int().unwrap(), 1);
        assert_eq!(Data::Bool(false).as_int().unwrap(), 0);
        assert_eq!(Data::Int(42).as_int().unwrap(), 42);
        assert_eq!(Data::Float(3.7).as_int().unwrap(), 3);
        assert!(Data::Bytes(vec![]).as_int().is_err());
    }

    #[test]
    fn test_data_len() {
        assert_eq!(Data::Bytes(b"test".to_vec()).len().unwrap(), 4);
        assert!(Data::Bool(true).len().is_err());
        assert!(Data::Int(42).len().is_err());
        assert!(Data::Float(2.71).len().is_err());
    }

    // -- Read methods --

    #[test]
    fn test_read_clear() {
        let origin = test_origin();
        let mut read = Read::new();
        read.add_fastq(1, b"read1", b"ACGT", b"IIII", origin, 0);
        assert!(!read.str_mappings.is_empty());
        read.clear();
        assert!(read.str_mappings.is_empty());
    }

    #[test]
    fn test_read_add_fastq_parts_with_qual() {
        let origin = test_origin();
        let mut read = Read::new();
        read.add_fastq_parts(1, Some(b"read1"), b"ACGT", Some(b"IIII"), origin, 0);
        assert!(read.str_mappings(StrType::Name(1)).is_some());
        let seq = read.str_mappings(StrType::Seq(1)).unwrap();
        assert_eq!(seq.string(), b"ACGT");
        assert_eq!(seq.qual(), Some(b"IIII".as_slice()));
    }

    #[test]
    fn test_read_add_fastq_parts_no_name() {
        let origin = test_origin();
        let mut read = Read::new();
        read.add_fastq_parts(1, None, b"ACGT", None, origin, 0);
        assert!(read.str_mappings(StrType::Name(1)).is_none());
        assert!(read.str_mappings(StrType::Seq(1)).is_some());
    }

    #[test]
    fn test_read_cut() {
        let origin = test_origin();
        let mut read = Read::new();
        read.add_fastq(1, b"read1", b"ACGTACGT", b"IIIIIIII", origin, 0);
        read.cut(
            StrType::Seq(1),
            InlineString::new(b"*"),
            Some(InlineString::new(b"left")),
            Some(InlineString::new(b"right")),
            4,
        )
        .unwrap();
        assert_eq!(
            read.substring(StrType::Seq(1), InlineString::new(b"left"))
                .unwrap(),
            b"ACGT"
        );
        assert_eq!(
            read.substring(StrType::Seq(1), InlineString::new(b"right"))
                .unwrap(),
            b"ACGT"
        );
    }

    #[test]
    fn test_read_set() {
        let origin = test_origin();
        let mut read = Read::new();
        read.add_fastq(1, b"read1", b"ACGTACGT", b"IIIIIIII", origin, 0);
        read.cut(
            StrType::Seq(1),
            InlineString::new(b"*"),
            Some(InlineString::new(b"left")),
            Some(InlineString::new(b"right")),
            4,
        )
        .unwrap();
        read.set(
            StrType::Seq(1),
            InlineString::new(b"left"),
            b"NN",
            Some(b"!!"),
        )
        .unwrap();
        assert_eq!(
            read.substring(StrType::Seq(1), InlineString::new(b"left"))
                .unwrap(),
            b"NN"
        );
    }

    #[test]
    fn test_read_trim() {
        let origin = test_origin();
        let mut read = Read::new();
        read.add_fastq(1, b"read1", b"ACGTACGT", b"IIIIIIII", origin, 0);
        read.cut(
            StrType::Seq(1),
            InlineString::new(b"*"),
            Some(InlineString::new(b"left")),
            Some(InlineString::new(b"right")),
            4,
        )
        .unwrap();
        read.trim(StrType::Seq(1), InlineString::new(b"left"))
            .unwrap();
        let seq = read.str_mappings(StrType::Seq(1)).unwrap();
        assert_eq!(seq.string(), b"ACGT");
    }

    #[test]
    fn test_read_remove_internal() {
        let origin = test_origin();
        let mut read = Read::new();
        read.add_fastq(1, b"read1", b"ACGTACGT", b"IIIIIIII", origin, 0);
        let sm = read.str_mappings_mut(StrType::Seq(1)).unwrap();
        sm.add_mapping(Some(InlineString::new(b"_internal")), 0, 4);
        sm.add_mapping(Some(InlineString::new(b"public")), 4, 4);
        read.remove_internal();
        let sm = read.str_mappings(StrType::Seq(1)).unwrap();
        assert!(sm.mapping(InlineString::new(b"_internal")).is_none());
        assert!(sm.mapping(InlineString::new(b"public")).is_some());
    }

    #[test]
    fn test_read_has_names() {
        let origin = test_origin();
        let mut read = Read::new();
        read.add_fastq(1, b"read1", b"ACGT", b"IIII", origin, 0);
        let label = crate::expr::Label::new(b"seq1.*").unwrap();
        assert!(read.has_names(&[crate::expr::LabelOrAttr::Label(label)]));
    }

    #[test]
    fn test_read_has_names_missing() {
        let read = Read::new();
        let label = crate::expr::Label::new(b"seq1.*").unwrap();
        assert!(!read.has_names(&[crate::expr::LabelOrAttr::Label(label)]));
    }

    #[test]
    fn test_read_data_mut() {
        let origin = test_origin();
        let mut read = Read::new();
        read.add_fastq(1, b"read1", b"ACGT", b"IIII", origin, 0);
        *read
            .data_mut(
                StrType::Seq(1),
                InlineString::new(b"*"),
                InlineString::new(b"score"),
            )
            .unwrap() = Data::Int(100);
        let val = read
            .data(
                StrType::Seq(1),
                InlineString::new(b"*"),
                InlineString::new(b"score"),
            )
            .unwrap();
        assert_eq!(val, &Data::Int(100));
    }

    #[test]
    fn test_read_to_fastq() {
        let origin = test_origin();
        let mut read = Read::new();
        read.add_fastq(1, b"read1", b"ACGT", b"IIII", origin, 0);
        let (name, seq, qual) = read.to_fastq(1).unwrap();
        assert_eq!(name, b"read1");
        assert_eq!(seq, b"ACGT");
        assert_eq!(qual, b"IIII");
    }

    #[test]
    fn test_read_to_json() {
        let origin = test_origin();
        let mut read = Read::new();
        read.add_fastq(1, b"read1", b"ACGT", b"IIII", origin, 0);
        let json = read.to_json();
        assert!(json.contains("ACGT"));
    }

    #[test]
    fn test_read_mapping() {
        let origin = test_origin();
        let mut read = Read::new();
        read.add_fastq(1, b"read1", b"ACGT", b"IIII", origin, 0);
        let m = read
            .mapping(StrType::Seq(1), InlineString::new(b"*"))
            .unwrap();
        assert_eq!(m.start, 0);
        assert_eq!(m.len, 4);
    }

    #[test]
    fn test_read_mapping_missing_str_type() {
        let read = Read::new();
        let result = read.mapping(StrType::Seq(1), InlineString::new(b"*"));
        assert!(result.is_err());
    }

    #[test]
    fn test_read_substring() {
        let origin = test_origin();
        let mut read = Read::new();
        read.add_fastq(1, b"read1", b"ACGT", b"IIII", origin, 0);
        assert_eq!(
            read.substring(StrType::Seq(1), InlineString::new(b"*"))
                .unwrap(),
            b"ACGT"
        );
    }

    #[test]
    fn test_read_substring_qual() {
        let origin = test_origin();
        let mut read = Read::new();
        read.add_fastq(1, b"read1", b"ACGT", b"IIII", origin, 0);
        assert_eq!(
            read.substring_qual(StrType::Seq(1), InlineString::new(b"*"))
                .unwrap(),
            Some(b"IIII".as_slice())
        );
    }

    #[test]
    fn test_read_intersect() {
        let origin = test_origin();
        let mut read = Read::new();
        read.add_fastq(1, b"read1", b"ACGTACGT", b"IIIIIIII", origin, 0);
        let sm = read.str_mappings_mut(StrType::Seq(1)).unwrap();
        sm.add_mapping(Some(InlineString::new(b"a")), 0, 5);
        sm.add_mapping(Some(InlineString::new(b"b")), 3, 5);
        read.intersect(
            StrType::Seq(1),
            InlineString::new(b"a"),
            InlineString::new(b"b"),
            Some(InlineString::new(b"inter")),
        )
        .unwrap();
        assert_eq!(
            read.substring(StrType::Seq(1), InlineString::new(b"inter"))
                .unwrap(),
            b"TA"
        );
    }

    #[test]
    fn test_read_union() {
        let origin = test_origin();
        let mut read = Read::new();
        read.add_fastq(1, b"read1", b"ACGTACGT", b"IIIIIIII", origin, 0);
        let sm = read.str_mappings_mut(StrType::Seq(1)).unwrap();
        sm.add_mapping(Some(InlineString::new(b"a")), 0, 3);
        sm.add_mapping(Some(InlineString::new(b"b")), 5, 3);
        read.union(
            StrType::Seq(1),
            InlineString::new(b"a"),
            InlineString::new(b"b"),
            Some(InlineString::new(b"uni")),
        )
        .unwrap();
        assert_eq!(
            read.substring(StrType::Seq(1), InlineString::new(b"uni"))
                .unwrap(),
            b"ACGTACGT"
        );
    }

    #[test]
    fn test_read_first_idx() {
        let origin = test_origin();
        let mut read = Read::new();
        read.add_fastq(1, b"read1", b"ACGT", b"IIII", origin, 5);
        assert_eq!(read.first_idx(), 5);
    }

    #[test]
    fn test_read_add_fastq_recycled() {
        let origin = test_origin();
        let mut read = Read::new();
        read.add_fastq(1, b"read1", b"ACGT", b"IIII", Arc::clone(&origin), 0);
        read.add_fastq_recycled(1, b"read2", b"TGCA", b"!!!!", origin, 1);
        let (name, seq, qual) = read.to_fastq(1).unwrap();
        assert_eq!(name, b"read2");
        assert_eq!(seq, b"TGCA");
        assert_eq!(qual, b"!!!!");
    }

    // -- StrMappings set with size changes --

    #[test]
    fn test_str_mappings_set_same_length() {
        let origin = test_origin();
        let mut sm = StrMappings::new(b"ACGTACGT".to_vec(), origin, 0);
        sm.cut(
            InlineString::new(b"*"),
            Some(InlineString::new(b"left")),
            Some(InlineString::new(b"right")),
            4,
        )
        .unwrap();
        sm.set(InlineString::new(b"left"), b"NNNN", None).unwrap();
        assert_eq!(sm.string(), b"NNNNACGT");
    }

    #[test]
    fn test_str_mappings_set_shorter() {
        let origin = test_origin();
        let mut sm = StrMappings::new(b"ACGTACGT".to_vec(), origin, 0);
        sm.cut(
            InlineString::new(b"*"),
            Some(InlineString::new(b"left")),
            Some(InlineString::new(b"right")),
            4,
        )
        .unwrap();
        sm.set(InlineString::new(b"left"), b"NN", None).unwrap();
        assert_eq!(sm.string(), b"NNACGT");
    }

    #[test]
    fn test_str_mappings_set_longer() {
        let origin = test_origin();
        let mut sm = StrMappings::new(b"ACGTACGT".to_vec(), origin, 0);
        sm.cut(
            InlineString::new(b"*"),
            Some(InlineString::new(b"left")),
            Some(InlineString::new(b"right")),
            4,
        )
        .unwrap();
        sm.set(InlineString::new(b"left"), b"NNNNNN", None).unwrap();
        assert_eq!(sm.string(), b"NNNNNNACGT");
    }

    // -- StrMappings display --

    #[test]
    fn test_str_mappings_display() {
        let origin = test_origin();
        let sm = StrMappings::new(b"ACGT".to_vec(), origin, 0);
        let display = format!("{}", sm);
        assert!(display.contains("ACGT"));
    }

    // -- Intersection edge cases with same start/end --

    #[test]
    fn test_intersection_same_start_a_longer() {
        let m1 = Mapping::new(InlineString::new(b"a"), 0, 8);
        let m2 = Mapping::new(InlineString::new(b"b"), 0, 4);
        let inter = m1.intersect(&m2);
        assert!(matches!(inter, Intersection::BAOverlap(_)));
    }

    #[test]
    fn test_intersection_same_start_b_longer() {
        let m1 = Mapping::new(InlineString::new(b"a"), 0, 4);
        let m2 = Mapping::new(InlineString::new(b"b"), 0, 8);
        let inter = m1.intersect(&m2);
        assert!(matches!(inter, Intersection::ABOverlap(_)));
    }

    #[test]
    fn test_intersection_same_end_a_starts_first() {
        // a=[0,8), b=[4,8) -> same end, a_start < b_start -> ABOverlap
        let m1 = Mapping::new(InlineString::new(b"a"), 0, 8);
        let m2 = Mapping::new(InlineString::new(b"b"), 4, 4);
        let inter = m1.intersect(&m2);
        assert!(matches!(inter, Intersection::ABOverlap(_)));
    }

    #[test]
    fn test_intersection_same_end_b_starts_first() {
        // a=[4,8), b=[0,8) -> same end, a_start > b_start -> BAOverlap
        let m1 = Mapping::new(InlineString::new(b"a"), 4, 4);
        let m2 = Mapping::new(InlineString::new(b"b"), 0, 8);
        let inter = m1.intersect(&m2);
        assert!(matches!(inter, Intersection::BAOverlap(_)));
    }

    // -- StrMappings recycle --

    #[test]
    fn test_str_mappings_recycle() {
        let origin = test_origin();
        let mut sm = StrMappings::new(b"ACGT".to_vec(), origin, 0);
        sm.add_mapping(Some(InlineString::new(b"test")), 0, 4);
        let origin2 = test_origin();
        sm.recycle(b"TGCA".to_vec(), origin2, 1);
        assert_eq!(sm.string(), b"TGCA");
        assert_eq!(sm.idx, 1);
        assert_eq!(sm.mappings.len(), 1);
    }

    #[test]
    fn test_str_mappings_recycle_with_qual() {
        let origin = test_origin();
        let mut sm = StrMappings::new_with_qual(b"ACGT".to_vec(), b"IIII".to_vec(), origin, 0);
        let origin2 = test_origin();
        sm.recycle_with_qual(b"TGCA".to_vec(), b"!!!!".to_vec(), origin2, 1);
        assert_eq!(sm.string(), b"TGCA");
        assert_eq!(sm.qual(), Some(b"!!!!".as_slice()));
    }
}
