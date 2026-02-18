use block_aligner::{cigar::*, scan_block::*, scores::*};

use rustc_hash::{FxHashMap, FxHashSet};

use memchr::memmem;

use thread_local::*;

use std::cell::RefCell;
use std::marker::Send;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::RwLock;

use crate::graph::*;
use crate::inline_string::InlineString;
use crate::seed_search::*;
use crate::Patterns;

/// Pre-computed lookup table for fast Hamming matching.
/// Stores substitution IDs as InlineString (Copy, stack-allocated) to avoid allocation.
struct HammingLookup {
    /// Maps encoded sequence -> substitution_id as InlineString (Copy type)
    table: FxHashMap<u64, InlineString>,
    /// Pattern length (all patterns must be same length)
    pattern_len: usize,
}

impl HammingLookup {
    const NUCLEOTIDES: [u8; 4] = [b'A', b'C', b'G', b'T'];

    /// Encode a sequence as u64 (up to 8 bytes)
    #[inline]
    fn encode(seq: &[u8]) -> u64 {
        let mut key = 0u64;
        for (i, &b) in seq.iter().enumerate() {
            key |= (b as u64) << (i * 8);
        }
        key
    }

    /// Build lookup table with all mismatch variants up to max_mismatches.
    /// sub_ids contains the substitution ID for each pattern index.
    fn new<'a>(
        patterns: impl Iterator<Item = (usize, &'a [u8])>,
        sub_ids: &[InlineString],
        pattern_len: usize,
        max_mismatches: usize,
    ) -> Self {
        let mut table = FxHashMap::default();

        for (pattern_idx, pattern) in patterns {
            let sub_id = sub_ids
                .get(pattern_idx)
                .copied()
                .unwrap_or_else(|| InlineString::new(b""));

            // Add exact match
            table.entry(Self::encode(pattern)).or_insert(sub_id);

            // Add 1-mismatch variants
            if max_mismatches >= 1 {
                for i in 0..pattern.len() {
                    for &nuc in &Self::NUCLEOTIDES {
                        if nuc != pattern[i] {
                            let mut variant = pattern.to_vec();
                            variant[i] = nuc;
                            table.entry(Self::encode(&variant)).or_insert(sub_id);
                        }
                    }
                }
            }

            // Add 2-mismatch variants
            if max_mismatches >= 2 {
                for i in 0..pattern.len() {
                    for j in (i + 1)..pattern.len() {
                        for &nuc1 in &Self::NUCLEOTIDES {
                            if nuc1 != pattern[i] {
                                for &nuc2 in &Self::NUCLEOTIDES {
                                    if nuc2 != pattern[j] {
                                        let mut variant = pattern.to_vec();
                                        variant[i] = nuc1;
                                        variant[j] = nuc2;
                                        table.entry(Self::encode(&variant)).or_insert(sub_id);
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }

        Self { table, pattern_len }
    }

    /// Lookup a sequence, returns the substitution ID if found.
    #[inline]
    fn lookup(&self, seq: &[u8]) -> Option<InlineString> {
        if seq.len() != self.pattern_len {
            return None;
        }
        self.table.get(&Self::encode(seq)).copied()
    }
}

pub struct MatchAnyOp {
    required_names: Vec<LabelOrAttr>,
    label: Label,
    new_labels: [Option<Label>; 3],
    patterns: Patterns,
    max_literal_len: usize,
    all_literals: bool,
    match_type: MatchType,
    aligner: ThreadLocal<Option<RefCell<Box<dyn Aligner + Send>>>>,
    seed_searcher: Option<SeedSearchers>,
    /// Fast hash-based lookup for Hamming matching (when applicable)
    hamming_lookup: Option<HammingLookup>,
    // Per-thread match-distance histograms; each thread stores counts by
    // exact edit distance, and we aggregate across threads when queried.
    distance_counts: ThreadLocal<RwLock<Vec<usize>>>,
    // Total number of reads that reached this node (across all threads).
    total_attempts: AtomicUsize,
}

impl MatchAnyOp {
    const NAME: &'static str = "MatchAnyOp";

    /// Match any one of multiple patterns in an interval.
    ///
    /// Patterns can be arbitrary expressions, so you can use any existing labeled intervals or
    /// attributes as patterns.
    ///
    /// You can also include arbitrary extra attributes for each pattern. The corresponding attributes
    /// for the matched pattern will be stored into the input labeled interval.
    ///
    /// The transform expression must have one input label and the number of output labels is
    /// determined by the [`MatchType`].
    ///
    /// Example `transform_expr` for local-alignment-based pattern matching:
    /// `tr!(seq1.* -> seq1.before, seq1.aligned, seq1.after)`.
    /// The input labeled interval will get a new attribute (`seq1.*.my_patterns`) that is set to the pattern
    /// that is matched. If no pattern matches, then it will be set to false.
    pub fn new(transform_expr: TransformExpr, patterns: Patterns, match_type: MatchType) -> Self {
        let mut new_labels = [None, None, None];

        transform_expr.check_size(1, match_type.num_mappings(), Self::NAME);
        for (i, label) in new_labels
            .iter_mut()
            .take(match_type.num_mappings())
            .enumerate()
        {
            *label = transform_expr.after_label(i, Self::NAME);
        }
        transform_expr.check_same_str_type(Self::NAME);

        let seed_searcher = Self::get_searcher(&patterns, &match_type);
        let max_literal_len = patterns
            .iter_literals()
            .map(|(_, p)| p.len())
            .max()
            .unwrap_or(0);
        let min_literal_len = patterns
            .iter_literals()
            .map(|(_, p)| p.len())
            .min()
            .unwrap_or(0);
        let all_literals = patterns.iter_exprs().count() == 0;
        let mut required_names = vec![transform_expr.before(0).into()];
        required_names.extend(
            patterns
                .iter_exprs()
                .flat_map(|(_, e)| e.required_names().into_iter()),
        );

        // Build fast hash-based lookup for Hamming matching when:
        // 1. All patterns are literals of the same length
        // 2. Match type is Hamming with small mismatch count (≤2)
        // 3. Pattern count is reasonable (≤1000)
        let hamming_lookup = if let MatchType::Hamming(threshold) = match_type {
            let max_mismatches = max_literal_len.saturating_sub(threshold.get(max_literal_len));
            let pattern_count = patterns.iter_literals().count();

            if all_literals
                && max_literal_len == min_literal_len
                && max_mismatches <= 2
                && pattern_count <= 1000
                && max_literal_len > 0
            {
                // Extract substitution IDs from pattern attributes as InlineStrings
                let sub_ids: Vec<InlineString> = patterns
                    .patterns()
                    .iter()
                    .map(|p| {
                        p.attrs()
                            .iter()
                            .find(|d| matches!(d, Data::Bytes(_)))
                            .and_then(|d| {
                                if let Data::Bytes(b) = d {
                                    // Convert to InlineString (up to 24 bytes)
                                    if b.len() <= 24 {
                                        Some(InlineString::new(b))
                                    } else {
                                        None
                                    }
                                } else {
                                    None
                                }
                            })
                            .unwrap_or_else(|| InlineString::new(b""))
                    })
                    .collect();

                Some(HammingLookup::new(
                    patterns.iter_literals(),
                    &sub_ids,
                    max_literal_len,
                    max_mismatches,
                ))
            } else {
                None
            }
        } else {
            None
        };

        Self {
            required_names,
            label: transform_expr.before(0),
            new_labels,
            patterns,
            max_literal_len,
            all_literals,
            match_type,
            aligner: ThreadLocal::new(),
            seed_searcher,
            hamming_lookup,
            distance_counts: ThreadLocal::new(),
            total_attempts: AtomicUsize::new(0),
        }
    }

    fn get_searcher(patterns: &Patterns, match_type: &MatchType) -> Option<SeedSearchers> {
        // For Edit/Hamming match types with small pattern sets, exhaustive search
        // (e.g. Myers bit-vector) is faster than building/querying the k-mer index.
        // Only skip seeding for these approximate match types; Exact and alignment
        // types should always use seeding when available.
        const MIN_PATTERNS_FOR_SEEDING: usize = 4;
        if matches!(
            match_type,
            MatchType::Edit(_)
                | MatchType::EditPrefix(_)
                | MatchType::EditSuffix(_)
                | MatchType::EditSearch(_)
                | MatchType::EditBoundedMatch { .. }
                | MatchType::Hamming(_)
                | MatchType::HammingPrefix(_)
                | MatchType::HammingSuffix(_)
                | MatchType::HammingSearch(_)
                | MatchType::HammingBoundedMatch { .. }
        ) && patterns.iter_literals().count() < MIN_PATTERNS_FOR_SEEDING
        {
            return None;
        }

        let min_len = patterns
            .iter_literals()
            .map(|(_, p)| p.len())
            .min()
            .unwrap_or(0);
        let k = match_type.k(min_len);

        use SeedSearchers::*;
        let res = match k {
            0..=1 => return None,
            2 => SmallSearcher::<2>::new(patterns.iter_literals()).map(Small2),
            3 => SmallSearcher::<3>::new(patterns.iter_literals()).map(Small3),
            4 => SmallSearcher::<4>::new(patterns.iter_literals()).map(Small4),
            5 => SmallSearcher::<5>::new(patterns.iter_literals()).map(Small5),
            6 => SmallSearcher::<6>::new(patterns.iter_literals()).map(Small6),
            _ => Err(()),
        };

        if let Ok(s) = res {
            Some(s)
        } else {
            Some(General(GeneralSearcher::new(patterns.iter_literals(), k)))
        }
    }

    #[inline]
    fn record_distance(&self, pattern_len: usize, matches: usize) {
        if pattern_len == 0 {
            return;
        }
        let distance = pattern_len.saturating_sub(matches);
        let cell = self.distance_counts.get_or(|| RwLock::new(Vec::new()));
        let mut counts = cell.write().unwrap();
        if distance >= counts.len() {
            counts.resize(distance + 1, 0);
        }
        counts[distance] += 1;
    }

    /// Human-readable label for statistics, e.g. "seq1.brc".
    pub fn stats_label(&self) -> String {
        format!("{}.{}", self.label.str_type, self.label.label)
    }
}

impl<T: crate::trace::Trace> GraphNode<T> for MatchAnyOp {
    fn run_inner(&self, mut reads: Vec<Read>) -> Result<(Option<Vec<Read>>, bool)> {
        // Count how many reads reach this node in this batch.
        self.total_attempts
            .fetch_add(reads.len(), Ordering::Relaxed);

        // Access thread-local aligner once per batch
        use MatchType::*;
        let aligner_cell = self.aligner.get_or(|| {
            let init_len = if self.max_literal_len > 0 {
                self.max_literal_len * 2
            } else {
                // Heuristic since we don't know text length yet, use reasonable default
                512
            };

            match self.match_type {
                GlobalAln(_) => Some(RefCell::new(Box::new(GlobalLocalAligner::<false>::new(
                    init_len,
                )))),
                LocalAln { .. } => Some(RefCell::new(Box::new(GlobalLocalAligner::<true>::new(
                    init_len,
                )))),
                PrefixAln { .. } => Some(RefCell::new(Box::new(PrefixSuffixAligner::<true>::new(
                    init_len,
                )))),
                SuffixAln { .. } => Some(RefCell::new(Box::new(
                    PrefixSuffixAligner::<false>::new(init_len),
                ))),
                _ => None,
            }
        });

        let additional = |identity: f64, pattern_len: usize| {
            ((1.0 - identity).max(0.0) * (pattern_len as f64)).ceil() as usize
        };

        for read in &mut reads {
            let text = read
                .substring(self.label.str_type, self.label.label)
                .map_err(|e| Error::NameError {
                    source: e,
                    read: read.clone(),
                    context: Self::NAME,
                })?;

            // Fast path: use pre-computed hash lookup for Hamming matching
            if let Some(ref lookup) = self.hamming_lookup {
                let pattern_len = lookup.pattern_len;
                match lookup.lookup(text) {
                    Some(sub_id) => {
                        // Fast path matched
                        let mapping = read
                            .mapping_mut(self.label.str_type, self.label.label)
                            .unwrap();

                        // Set sub and ambig attributes
                        let sub_bytes: Vec<u8> = sub_id.bytes().collect();
                        *mapping.data_mut(InlineString::new(b"sub")) = Data::Bytes(sub_bytes);
                        *mapping.data_mut(InlineString::new(b"ambig")) =
                            Data::Bytes(b"false".to_vec());

                        // For Hamming match, num_mappings() is 1
                        let start = mapping.start;
                        let str_mappings = read.str_mappings_mut(self.label.str_type).unwrap();
                        str_mappings.add_mapping(
                            self.new_labels[0].as_ref().map(|l| l.label),
                            start,
                            pattern_len,
                        );
                        continue; // Skip slow path
                    }
                    None => {
                        // Fast path no match - set up no-match result
                        let (start, len) = {
                            let mapping =
                                read.mapping(self.label.str_type, self.label.label).unwrap();
                            (mapping.start, mapping.len)
                        };

                        if let Some(new_label) = &self.new_labels[0] {
                            let str_mappings = read.str_mappings_mut(self.label.str_type).unwrap();
                            str_mappings.add_mapping(Some(new_label.label), start, len);
                        }

                        let mapping = read
                            .mapping_mut(self.label.str_type, self.label.label)
                            .unwrap();

                        // Set empty/ambig attributes for no-match
                        *mapping.data_mut(InlineString::new(b"sub")) = Data::Bytes(Vec::new());
                        *mapping.data_mut(InlineString::new(b"ambig")) =
                            Data::Bytes(b"true".to_vec());
                        continue; // Skip slow path
                    }
                }
            }

            let mut seed_hits = FxHashSet::default();

            if let Some(seed_searcher) = &self.seed_searcher {
                let (text_slice, text_offset, use_i) = match self.match_type {
                    Exact => (text, 0, false),
                    ExactPrefix => (&text[..text.len().min(self.max_literal_len)], 0, false),
                    ExactSuffix => {
                        let offset = text.len().saturating_sub(self.max_literal_len);
                        (&text[offset..], offset, false)
                    }
                    ExactSearch => (text, 0, true),
                    ExactBoundedMatch { from, to } => {
                        let to = text.len().min(to);
                        (&text[from..to], 0, false)
                    }
                    Hamming(_) => (text, 0, false),
                    HammingPrefix(_) => (&text[..text.len().min(self.max_literal_len)], 0, false),
                    HammingSuffix(_) => {
                        let offset = text.len().saturating_sub(self.max_literal_len);
                        (&text[offset..], offset, false)
                    }
                    HammingSearch(_) => (text, 0, true),
                    HammingBoundedMatch {
                        threshold: _,
                        from,
                        to,
                    } => {
                        let to = text.len().min(to);
                        (&text[from..to], 0, false)
                    }
                    GlobalAln(_) => (text, 0, false),
                    LocalAln { .. } => (text, 0, true),
                    PrefixAln { identity, .. } => (
                        &text[..text.len().min(
                            self.max_literal_len + additional(identity, self.max_literal_len),
                        )],
                        0,
                        false,
                    ),
                    SuffixAln { identity, .. } => {
                        let offset = text.len().saturating_sub(
                            self.max_literal_len + additional(identity, self.max_literal_len),
                        );
                        (&text[offset..], offset, false)
                    }
                    Edit(_) => (text, 0, false),
                    EditPrefix(t) => {
                        let max_edits = t.get(self.max_literal_len);
                        (
                            &text[..text.len().min(self.max_literal_len + max_edits)],
                            0,
                            false,
                        )
                    }
                    EditSuffix(t) => {
                        let max_edits = t.get(self.max_literal_len);
                        let offset = text.len().saturating_sub(self.max_literal_len + max_edits);
                        (&text[offset..], offset, false)
                    }
                    EditSearch(_) => (text, 0, true),
                    EditBoundedMatch {
                        threshold: _,
                        from,
                        to,
                    } => {
                        let to = text.len().min(to);
                        (&text[from..to], 0, false)
                    }
                };

                seed_searcher.search(
                    text_slice,
                    |SeedMatch {
                         pattern_idx,
                         pattern_i,
                         text_i,
                     }| {
                        let text_i = if use_i {
                            Some(((text_offset + text_i) as isize) - (pattern_i as isize))
                        } else {
                            None
                        };
                        seed_hits.insert((pattern_idx, text_i));
                    },
                );
            } else {
                seed_hits.extend(self.patterns.iter_literals().map(|(i, _)| (i, None)));
            }

            if !self.all_literals {
                seed_hits.extend(self.patterns.iter_exprs().map(|(i, _)| (i, None)));
            }

            let mut max_matches = 0;
            let mut max_pattern_len = 0;
            let mut max_pattern = None;
            let mut max_pattern_idx = usize::MAX;
            let mut max_cut_pos1 = 0;
            let mut max_cut_pos2 = 0;
            let mut multimatches = false;

            for (pattern_idx, text_i) in seed_hits {
                let pattern = &self.patterns.patterns()[pattern_idx];
                let pattern_str_cow = pattern.get(read).map_err(|e| Error::NameError {
                    source: e,
                    read: read.clone(),
                    context: Self::NAME,
                })?;
                let pattern_str: &[u8] = &pattern_str_cow;
                let pattern_len = pattern_str.len();

                if max_matches > pattern_len {
                    continue;
                }

                let matches = match self.match_type {
                    Exact => {
                        if text == pattern_str {
                            Some((pattern_len, pattern_len, 0))
                        } else {
                            None
                        }
                    }
                    ExactPrefix => {
                        if pattern_len <= text.len() && &text[..pattern_len] == pattern_str {
                            Some((pattern_len, pattern_len, 0))
                        } else {
                            None
                        }
                    }
                    ExactSuffix => {
                        if pattern_len <= text.len()
                            && &text[text.len() - pattern_len..] == pattern_str
                        {
                            Some((pattern_len, text.len() - pattern_len, 0))
                        } else {
                            None
                        }
                    }
                    ExactSearch => {
                        let (text_start, text_end) = if let Some(text_i) = text_i {
                            (
                                text_i.max(0) as usize,
                                text.len().min((text_i + (pattern_len as isize)) as usize),
                            )
                        } else {
                            (0, text.len())
                        };
                        let text_around = &text[text_start..text_end];
                        memmem::find(text_around, pattern_str)
                            .map(|i| (pattern_len, text_start + i, text_start + i + pattern_len))
                    }
                    ExactBoundedMatch { from, to } => {
                        let to = text.len().min(to);
                        let text_around = &text[from..=to];
                        memmem::find(text_around, pattern_str)
                            .map(|i| (pattern_len, from + i, from + i + pattern_len))
                    }
                    Hamming(t) => {
                        let t = t.get(pattern_len);
                        hamming(text, pattern_str, t).map(|m| (m, pattern_len, 0))
                    }
                    HammingPrefix(t) => {
                        if pattern_len <= text.len() {
                            let t = t.get(pattern_len);
                            hamming(&text[..pattern_len], pattern_str, t)
                                .map(|m| (m, pattern_len, 0))
                        } else {
                            None
                        }
                    }
                    HammingSuffix(t) => {
                        if pattern_len <= text.len() {
                            let t = t.get(pattern_len);
                            hamming(&text[text.len() - pattern_len..], pattern_str, t)
                                .map(|m| (m, text.len() - pattern_len, 0))
                        } else {
                            None
                        }
                    }
                    HammingSearch(t) => {
                        let t = t.get(pattern_len);
                        if let Some(text_i) = text_i {
                            // Seed hit gives us the exact position - just check that position
                            let text_start = text_i.max(0) as usize;
                            let text_end = text.len().min(text_start + pattern_len);
                            if text_end - text_start == pattern_len {
                                let text_slice = &text[text_start..text_end];
                                hamming(text_slice, pattern_str, t)
                                    .map(|m| (m, text_start, text_end))
                            } else {
                                None
                            }
                        } else {
                            // No seed hit - fall back to full search
                            hamming_search(text, pattern_str, t)
                        }
                    }
                    HammingBoundedMatch {
                        threshold: t,
                        from,
                        to,
                    } => {
                        let t = t.get(pattern_len);
                        // Use exclusive range - to is the max position, so we need to+1 for the slice
                        // but capped at text.len()
                        let to_exclusive = text.len().min(to + 1);
                        let text_around = &text[from..to_exclusive];
                        hamming_search(text_around, pattern_str, t)
                            .map(|(m, start_idx, end_idx)| (m, from + start_idx, from + end_idx))
                    }
                    GlobalAln(identity) => aligner_cell
                        .as_ref()
                        .unwrap()
                        .borrow_mut()
                        .align(text, pattern_str, identity, identity)
                        .map(|(m, _, end_idx)| (m, end_idx, 0)),
                    LocalAln { identity, overlap } => {
                        let a = additional(identity, pattern_len) as isize;
                        let (text_start, text_end) = if let Some(text_i) = text_i {
                            (
                                (text_i - a).max(0) as usize,
                                text.len()
                                    .min((text_i + (pattern_len as isize) + a) as usize),
                            )
                        } else {
                            (0, text.len())
                        };
                        let text_around = &text[text_start..text_end];
                        aligner_cell
                            .as_ref()
                            .unwrap()
                            .borrow_mut()
                            .align(text_around, pattern_str, identity, overlap)
                            .map(|(m, start_idx, end_idx)| {
                                (m, text_start + start_idx, text_start + end_idx)
                            })
                    }
                    PrefixAln { identity, overlap } => {
                        let a = additional(identity, pattern_len);
                        aligner_cell
                            .as_ref()
                            .unwrap()
                            .borrow_mut()
                            .align(
                                &text[..text.len().min(pattern_len + a)],
                                pattern_str,
                                identity,
                                overlap,
                            )
                            .map(|(m, _, end_idx)| (m, end_idx, 0))
                    }
                    SuffixAln { identity, overlap } => {
                        let a = additional(identity, pattern_len);
                        let text_start = text.len().saturating_sub(pattern_len + a);
                        aligner_cell
                            .as_ref()
                            .unwrap()
                            .borrow_mut()
                            .align(&text[text_start..], pattern_str, identity, overlap)
                            .map(|(m, start_idx, _)| (m, text_start + start_idx, 0))
                    }
                    Edit(t) => {
                        let max_edits = t.get(pattern_len);
                        edit_distance(text, pattern_str, max_edits).map(|m| (m, pattern_len, 0))
                    }
                    EditPrefix(t) => {
                        let max_edits = t.get(pattern_len);
                        edit_prefix(text, pattern_str, max_edits)
                            .map(|(m, end_pos)| (m, end_pos, 0))
                    }
                    EditSuffix(t) => {
                        let max_edits = t.get(pattern_len);
                        edit_suffix(text, pattern_str, max_edits)
                            .map(|(m, start_pos)| (m, start_pos, 0))
                    }
                    EditSearch(t) => {
                        let max_edits = t.get(pattern_len);
                        if let Some(text_i) = text_i {
                            // Seed hit - search around the seed position
                            let text_start = (text_i - (max_edits as isize)).max(0) as usize;
                            let text_end = text.len().min(
                                (text_i + (pattern_len as isize) + (max_edits as isize)) as usize,
                            );
                            if text_end > text_start {
                                let text_slice = &text[text_start..text_end];
                                edit_search(text_slice, pattern_str, max_edits).map(
                                    |(m, start_idx, end_idx)| {
                                        (m, text_start + start_idx, text_start + end_idx)
                                    },
                                )
                            } else {
                                None
                            }
                        } else {
                            // No seed hit - full search
                            edit_search(text, pattern_str, max_edits)
                        }
                    }
                    EditBoundedMatch {
                        threshold: t,
                        from,
                        to,
                    } => {
                        let max_edits = t.get(pattern_len);
                        let to_exclusive = text.len().min(to + 1);
                        let text_around = &text[from..to_exclusive];
                        edit_search(text_around, pattern_str, max_edits)
                            .map(|(m, start_idx, end_idx)| (m, from + start_idx, from + end_idx))
                    }
                };

                if let Some((matches, cut_pos1, cut_pos2)) = matches {
                    if matches > max_matches {
                        max_matches = matches;
                        max_pattern_len = pattern_len;
                        max_pattern = Some((pattern_str_cow, pattern.attrs()));
                        max_pattern_idx = pattern_idx;
                        max_cut_pos1 = cut_pos1;
                        max_cut_pos2 = cut_pos2;
                        multimatches = false;
                    } else if matches == max_matches && pattern_idx != max_pattern_idx {
                        multimatches = true;
                    }
                }
            }

            if let Some((pattern_str, pattern_attrs)) = max_pattern {
                self.record_distance(max_pattern_len, max_matches);
                let pattern_str = pattern_str.into_owned();
                let mapping = read
                    .mapping_mut(self.label.str_type, self.label.label)
                    .unwrap();

                if let Some(pattern_name) = self.patterns.pattern_name() {
                    *mapping.data_mut(pattern_name) = Data::Bytes(pattern_str);
                }

                if let Some(multimatch_name) = self.patterns.multimatch_name() {
                    let val = if multimatches {
                        b"true".to_vec()
                    } else {
                        b"false".to_vec()
                    };
                    *mapping.data_mut(multimatch_name) = Data::Bytes(val);
                }

                for (&attr, data) in self.patterns.attr_names().iter().zip(pattern_attrs) {
                    *mapping.data_mut(attr) = data.clone();
                }

                match self.match_type.num_mappings() {
                    1 => {
                        let start = mapping.start;
                        let str_mappings = read.str_mappings_mut(self.label.str_type).unwrap();
                        str_mappings.add_mapping(
                            self.new_labels[0].as_ref().map(|l| l.label),
                            start,
                            max_cut_pos1,
                        );
                    }
                    2 => {
                        read.cut(
                            self.label.str_type,
                            self.label.label,
                            self.new_labels[0].as_ref().map(|l| l.label),
                            self.new_labels[1].as_ref().map(|l| l.label),
                            max_cut_pos1 as isize,
                        )
                        .unwrap_or_else(|e| panic!("Error in {}: {e}", Self::NAME));
                    }
                    3 => {
                        let offset = mapping.start;
                        let mapping_len = mapping.len;

                        let str_mappings = read.str_mappings_mut(self.label.str_type).unwrap();
                        str_mappings.add_mapping(
                            self.new_labels[0].as_ref().map(|l| l.label),
                            offset,
                            max_cut_pos1,
                        );
                        str_mappings.add_mapping(
                            self.new_labels[1].as_ref().map(|l| l.label),
                            offset + max_cut_pos1,
                            max_cut_pos2 - max_cut_pos1,
                        );
                        str_mappings.add_mapping(
                            self.new_labels[2].as_ref().map(|l| l.label),
                            offset + max_cut_pos2,
                            mapping_len - max_cut_pos2,
                        );
                    }
                    _ => unreachable!(),
                }
            } else {
                let (start, len) = {
                    let mapping = read.mapping(self.label.str_type, self.label.label).unwrap();
                    (mapping.start, mapping.len)
                };

                // Pass-through the label if it's a 1-to-1 transform (e.g. MatchType::Exact)
                if self.match_type.num_mappings() == 1 {
                    if let Some(new_label) = &self.new_labels[0] {
                        let str_mappings = read.str_mappings_mut(self.label.str_type).unwrap();
                        str_mappings.add_mapping(Some(new_label.label), start, len);
                    }
                }

                let mapping = read
                    .mapping_mut(self.label.str_type, self.label.label)
                    .unwrap();

                // Reset pattern name (unused by seqproc map) on no-match
                if let Some(pattern_name) = self.patterns.pattern_name() {
                    *mapping.data_mut(pattern_name) = Data::Bytes(Vec::new());
                }

                // For seqproc's map(), `ambig` is used to derive the boolean `MAPPED = !ambig`.
                // On an unmatched read we want MAPPED == false so that:
                //   * the fallback graph (e.g. pad_to) runs, and
                //   * the mapping graph that dereferences `.sub` is NOT executed.
                // Since `expect_bool` treats non-empty, non-"false" bytes as true,
                // we store "true" here so that `!ambig` evaluates to false.
                if let Some(multimatch_name) = self.patterns.multimatch_name() {
                    *mapping.data_mut(multimatch_name) = Data::Bytes(b"true".to_vec());
                }

                // Initialize any pattern attributes; `sub` is left as empty bytes for no-match.
                for &attr in self.patterns.attr_names() {
                    let name = attr.as_str().as_bytes();
                    if name == b"sub" {
                        *mapping.data_mut(attr) = Data::Bytes(Vec::new());
                    } else if name == b"ambig" {
                        *mapping.data_mut(attr) = Data::Bytes(b"true".to_vec());
                    } else {
                        *mapping.data_mut(attr) = Data::Bytes(Vec::new());
                    }
                }

                // Force-create defaults for `sub`/`ambig` if they were not in attr_names.
                *mapping.data_mut(InlineString::new(b"sub")) = Data::Bytes(Vec::new());
                *mapping.data_mut(InlineString::new(b"ambig")) = Data::Bytes(b"true".to_vec());
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

    fn match_distance_counts(&self) -> Option<MatchDistanceCounts> {
        let mut totals: Vec<usize> = Vec::new();

        for local in self.distance_counts.iter() {
            let local = local.read().unwrap();
            if local.len() > totals.len() {
                totals.resize(local.len(), 0);
            }
            for (d, &count) in local.iter().enumerate() {
                totals[d] += count;
            }
        }

        // Trim trailing zeros to keep the internal representation compact.
        while totals.last().copied() == Some(0) {
            totals.pop();
        }

        Some(MatchDistanceCounts {
            label: self.stats_label(),
            counts: totals,
            total: self.total_attempts.load(Ordering::Relaxed),
        })
    }
}

fn hamming(a: &[u8], b: &[u8], threshold: usize) -> Option<usize> {
    if a.len() != b.len() {
        return None;
    }

    let a_ptr = a.as_ptr();
    let b_ptr = b.as_ptr();
    let n = a.len();
    let mut res = 0;
    let mut i = 0;

    unsafe {
        while i < (n / 8) * 8 {
            let a_word = std::ptr::read_unaligned(a_ptr.add(i) as *const u64);
            let b_word = std::ptr::read_unaligned(b_ptr.add(i) as *const u64);

            let xor = a_word ^ b_word;
            let or1 = xor | (xor >> 1);
            let or2 = or1 | (or1 >> 2);
            let or3 = or2 | (or2 >> 4);
            let mask = or3 & 0x0101010101010101u64;
            res += mask.count_ones() as usize;

            i += 8;
        }

        if i < n {
            let a_word = read_rest_u64(a_ptr.add(i), n - i);
            let b_word = read_rest_u64(b_ptr.add(i), n - i);

            let xor = a_word ^ b_word;
            let or1 = xor | (xor >> 1);
            let or2 = or1 | (or1 >> 2);
            let or3 = or2 | (or2 >> 4);
            let mask = or3 & 0x0101010101010101u64;
            res += mask.count_ones() as usize;
        }
    }

    let matches = n - res;

    if matches >= threshold {
        Some(matches)
    } else {
        None
    }
}

unsafe fn read_rest_u64(ptr: *const u8, len: usize) -> u64 {
    let addr = ptr as usize;
    let start_page = addr >> 12;
    let end_page = (addr + 7) >> 12;

    if start_page == end_page {
        std::ptr::read_unaligned(ptr as *const u64) & ((1u64 << (len * 8)) - 1)
    } else {
        let mut res = 0u64;
        let mut i = 0;

        while i < len {
            res |= (*ptr.add(i) as u64) << (i * 8);
            i += 1;
        }

        res
    }
}

fn hamming_search(a: &[u8], b: &[u8], threshold: usize) -> Option<(usize, usize, usize)> {
    let mut best_match = None;

    for (i, w) in a.windows(b.len()).enumerate() {
        if let Some(matches) = hamming(w, b, threshold) {
            if let Some((best_matches, _, _)) = best_match {
                if matches <= best_matches {
                    continue;
                }
            }

            best_match = Some((matches, i, i + b.len()));
        }
    }

    best_match
}

/// Compute edit distance (Levenshtein distance) between two sequences.
/// Returns the number of matches (len - edits) if within threshold, None otherwise.
/// Uses Myers' bit-vector algorithm for sequences up to 64bp, falls back to DP otherwise.
fn edit_distance(text: &[u8], pattern: &[u8], max_edits: usize) -> Option<usize> {
    let m = pattern.len();
    let n = text.len();

    if m == 0 {
        return if n <= max_edits { Some(m) } else { None };
    }
    if n == 0 {
        return if m <= max_edits { Some(0) } else { None };
    }

    // For full match, lengths should be similar within edit distance
    let len_diff = n.abs_diff(m);
    if len_diff > max_edits {
        return None;
    }

    // Use Myers' bit-vector for patterns up to 64bp
    if m <= 64 {
        edit_distance_myers(text, pattern, max_edits)
    } else {
        edit_distance_dp(text, pattern, max_edits)
    }
}

/// Myers' bit-vector algorithm for edit distance (patterns up to 64bp)
fn edit_distance_myers(text: &[u8], pattern: &[u8], max_edits: usize) -> Option<usize> {
    let m = pattern.len();
    let _n = text.len();

    // Build pattern bitmasks for each character
    let mut peq = [0u64; 256];
    for (i, &c) in pattern.iter().enumerate() {
        peq[c as usize] |= 1u64 << i;
    }

    // Initialize bit vectors
    let mut pv: u64 = !0u64; // all 1s
    let mut mv: u64 = 0u64; // all 0s
    let mut score = m;
    let high_bit = 1u64 << (m - 1);

    for &c in text {
        let eq = peq[c as usize];
        let xv = eq | mv;
        let xh = ((eq & pv).wrapping_add(pv)) ^ pv | eq;

        let ph = mv | !(xh | pv);
        let mh = pv & xh;

        // Update score
        if (ph & high_bit) != 0 {
            score += 1;
        }
        if (mh & high_bit) != 0 {
            score -= 1;
        }

        // Shift for next iteration
        pv = (mh << 1) | !(xv | (ph << 1));
        mv = (ph << 1) & xv;
    }

    if score <= max_edits {
        Some(m.saturating_sub(score))
    } else {
        None
    }
}

/// Standard DP algorithm for edit distance (for patterns > 64bp)
fn edit_distance_dp(text: &[u8], pattern: &[u8], max_edits: usize) -> Option<usize> {
    let m = pattern.len();
    let n = text.len();

    // Use two rows for space efficiency
    let mut prev = vec![0usize; m + 1];
    let mut curr = vec![0usize; m + 1];

    // Initialize first row
    for (j, val) in prev.iter_mut().enumerate() {
        *val = j;
    }

    for i in 1..=n {
        curr[0] = i;
        let mut min_in_row = curr[0];

        for j in 1..=m {
            let cost = if text[i - 1] == pattern[j - 1] { 0 } else { 1 };
            curr[j] = (prev[j - 1] + cost).min(prev[j] + 1).min(curr[j - 1] + 1);
            min_in_row = min_in_row.min(curr[j]);
        }

        // Early termination if minimum in row exceeds threshold
        if min_in_row > max_edits {
            return None;
        }

        std::mem::swap(&mut prev, &mut curr);
    }

    let edits = prev[m];
    if edits <= max_edits {
        Some(m.saturating_sub(edits))
    } else {
        None
    }
}

/// Search for the best edit distance match of pattern in text.
/// Returns (matches, start_idx, end_idx) for the best match within threshold.
/// Uses semi-global alignment where gaps at text boundaries are free.
fn edit_search(text: &[u8], pattern: &[u8], max_edits: usize) -> Option<(usize, usize, usize)> {
    let m = pattern.len();
    let n = text.len();

    if m == 0 || n == 0 {
        return None;
    }

    // Use Myers' bit-vector semi-global search for patterns up to 64bp
    if m <= 64 {
        edit_search_myers(text, pattern, max_edits)
    } else {
        edit_search_dp(text, pattern, max_edits)
    }
}

/// Myers' bit-vector algorithm for semi-global edit distance search (patterns up to 64bp).
/// Uses a forward pass to find the best end position, then a reverse DP pass to find
/// the exact start position.
fn edit_search_myers(
    text: &[u8],
    pattern: &[u8],
    max_edits: usize,
) -> Option<(usize, usize, usize)> {
    let m = pattern.len();

    // Build pattern bitmasks
    let mut peq = [0u64; 256];
    for (i, &c) in pattern.iter().enumerate() {
        peq[c as usize] |= 1u64 << i;
    }

    let mut pv: u64 = !0u64;
    let mut mv: u64 = 0u64;
    let mut score = m;
    let high_bit = 1u64 << (m - 1);

    let mut best_score = usize::MAX;
    let mut best_end = 0;

    for (i, &c) in text.iter().enumerate() {
        let eq = peq[c as usize];
        let xv = eq | mv;
        let xh = ((eq & pv).wrapping_add(pv)) ^ pv | eq;

        let ph = mv | !(xh | pv);
        let mh = pv & xh;

        if (ph & high_bit) != 0 {
            score += 1;
        }
        if (mh & high_bit) != 0 {
            score -= 1;
        }

        if score <= max_edits && score < best_score {
            best_score = score;
            best_end = i + 1;
        }

        pv = (mh << 1) | !(xv | (ph << 1));
        mv = (ph << 1) & xv;
    }

    if best_score > max_edits {
        return None;
    }

    // Reverse DP to find exact start position
    let start = find_start_reverse_dp(&text[..best_end], pattern, best_score);

    Some((m.saturating_sub(best_score), start, best_end))
}

/// Find the exact start position of a semi-global alignment by running a reverse DP.
/// Given that the best alignment ends at text[..end_pos] with `edits` edits,
/// align the reversed pattern against the reversed text suffix to find where
/// the alignment begins.
#[allow(clippy::needless_range_loop)]
fn find_start_reverse_dp(text: &[u8], pattern: &[u8], edits: usize) -> usize {
    let m = pattern.len();
    let n = text.len();

    // Reverse DP: align reversed pattern against reversed text (semi-global)
    let mut prev = vec![0usize; m + 1];
    let mut curr = vec![0usize; m + 1];

    for j in 0..=m {
        prev[j] = j;
    }

    let mut best_score = m;
    let mut best_start_from_end = 0;

    for i in 1..=n {
        curr[0] = 0; // Free gaps at start of reversed text (= free gaps at end of original)
        for j in 1..=m {
            let cost = if text[n - i] == pattern[m - j] { 0 } else { 1 };
            curr[j] = (prev[j - 1] + cost).min(prev[j] + 1).min(curr[j - 1] + 1);
        }
        if curr[m] <= edits && curr[m] <= best_score {
            best_score = curr[m];
            best_start_from_end = i;
        }
        std::mem::swap(&mut prev, &mut curr);
    }

    n - best_start_from_end
}

/// DP-based semi-global edit distance search for longer patterns.
/// Uses a forward pass to find the best end position, then a reverse DP pass to find
/// the exact start position.
#[allow(clippy::needless_range_loop)]
fn edit_search_dp(text: &[u8], pattern: &[u8], max_edits: usize) -> Option<(usize, usize, usize)> {
    let m = pattern.len();
    let n = text.len();

    // DP with free gaps at text boundaries (semi-global)
    let mut prev = vec![0usize; m + 1];
    let mut curr = vec![0usize; m + 1];

    // Initialize: gaps in pattern cost, gaps in text at start are free
    for j in 0..=m {
        prev[j] = j;
    }

    let mut best_score = usize::MAX;
    let mut best_end = 0;

    for i in 1..=n {
        curr[0] = 0; // Free gaps at text start

        for j in 1..=m {
            let cost = if text[i - 1] == pattern[j - 1] { 0 } else { 1 };
            curr[j] = (prev[j - 1] + cost).min(prev[j] + 1).min(curr[j - 1] + 1);
        }

        // Check if this is a valid end position (free gaps at text end)
        if curr[m] <= max_edits && curr[m] < best_score {
            best_score = curr[m];
            best_end = i;
        }

        std::mem::swap(&mut prev, &mut curr);
    }

    if best_score > max_edits {
        return None;
    }

    // Reverse DP to find exact start position
    let start = find_start_reverse_dp(&text[..best_end], pattern, best_score);

    Some((m.saturating_sub(best_score), start, best_end))
}

/// Edit distance for prefix matching - pattern should match a prefix of text.
/// Uses a DP where gaps at the text start DO cost (alignment must begin at position 0),
/// but gaps at the text end are free (the match can end anywhere).
/// Returns (matches, end_position) where matches = pattern_len - edits.
#[allow(clippy::needless_range_loop)]
fn edit_prefix(text: &[u8], pattern: &[u8], max_edits: usize) -> Option<(usize, usize)> {
    let m = pattern.len();

    if m == 0 {
        return Some((0, 0));
    }

    // Only need to scan up to m + max_edits text positions
    let search_end = (m + max_edits).min(text.len());
    let text_prefix = &text[..search_end];
    let n = text_prefix.len();

    let mut prev = vec![0usize; m + 1];
    let mut curr = vec![0usize; m + 1];

    // Initialize: dp[0][j] = j (cost to match first j pattern chars with empty text prefix)
    for j in 0..=m {
        prev[j] = j;
    }

    let mut best_score = prev[m]; // aligning empty text against full pattern = m deletions
    let mut best_end = 0;

    for i in 1..=n {
        curr[0] = i; // Gaps at text start DO cost (unlike semi-global search)

        for j in 1..=m {
            let cost = if text_prefix[i - 1] == pattern[j - 1] {
                0
            } else {
                1
            };
            curr[j] = (prev[j - 1] + cost).min(prev[j] + 1).min(curr[j - 1] + 1);
        }

        // Free gaps at text end: check if full pattern is matched at this text position
        if curr[m] <= max_edits && curr[m] <= best_score {
            best_score = curr[m];
            best_end = i;
        }

        std::mem::swap(&mut prev, &mut curr);
    }

    if best_score > max_edits {
        return None;
    }

    Some((m.saturating_sub(best_score), best_end))
}

/// Edit distance for suffix matching - pattern should match a suffix of text.
/// Uses a DP where gaps at the text end DO cost (alignment must end at the last position),
/// but gaps at the text start are free (the match can begin anywhere).
/// Returns (matches, start_position) where matches = pattern_len - edits.
#[allow(clippy::needless_range_loop)]
fn edit_suffix(text: &[u8], pattern: &[u8], max_edits: usize) -> Option<(usize, usize)> {
    let m = pattern.len();
    let n = text.len();

    if m == 0 {
        return Some((0, n));
    }

    // Only need to scan the last m + max_edits text positions
    let search_start = n.saturating_sub(m + max_edits);
    let text_suffix = &text[search_start..];
    let sn = text_suffix.len();

    // Reverse DP on text_suffix and pattern:
    // Align reversed pattern against reversed text_suffix.
    // Free gaps at the start of reversed text (= free gaps at the END of original text_suffix)
    // would be wrong -- we want suffix alignment where the match must reach the text end.
    // Instead: align reversed text_suffix against reversed pattern with:
    //   curr[0] = i (gaps at reversed-text start cost = gaps at original-text end cost)
    //   answer = min over i of dp[i][m] (free gaps at reversed-text end = free original-text start)
    // This is the mirror of edit_prefix.
    let mut prev = vec![0usize; m + 1];
    let mut curr = vec![0usize; m + 1];

    for j in 0..=m {
        prev[j] = j;
    }

    let mut best_score = prev[m];
    let mut best_start_from_end = 0;

    for i in 1..=sn {
        curr[0] = i; // Gaps at text end DO cost

        for j in 1..=m {
            // Traverse both text and pattern in reverse
            let cost = if text_suffix[sn - i] == pattern[m - j] {
                0
            } else {
                1
            };
            curr[j] = (prev[j - 1] + cost).min(prev[j] + 1).min(curr[j - 1] + 1);
        }

        // Free gaps at text start: check if full pattern is matched at this text position
        if curr[m] <= max_edits && curr[m] <= best_score {
            best_score = curr[m];
            best_start_from_end = i;
        }

        std::mem::swap(&mut prev, &mut curr);
    }

    if best_score > max_edits {
        return None;
    }

    let start = search_start + (sn - best_start_from_end);
    Some((m.saturating_sub(best_score), start))
}

trait Aligner {
    fn align(
        &mut self,
        read: &[u8],
        pattern: &[u8],
        identity_threshold: f64,
        overlap_threshold: f64,
    ) -> Option<(usize, usize, usize)>;
}

struct GlobalLocalAligner<const LOCAL: bool> {
    read_padded: PaddedBytes,
    pattern_padded: PaddedBytes,
    matrix: NucMatrix,
    // always store trace
    block: Block<true, LOCAL, LOCAL, false>,
    cigar: Cigar,
    len: usize,
}

impl<const LOCAL: bool> GlobalLocalAligner<LOCAL> {
    const MIN_SIZE: usize = 32;
    const MAX_SIZE: usize = 512;
    const GAPS: Gaps = Gaps {
        open: -2,
        extend: -1,
    };

    pub fn new(len: usize) -> Self {
        let read_padded = PaddedBytes::new::<NucMatrix>(len, Self::MAX_SIZE);
        let pattern_padded = PaddedBytes::new::<NucMatrix>(len, Self::MAX_SIZE);
        let matrix = NucMatrix::new_simple(1, -1);

        let block = Block::<true, LOCAL, LOCAL, false>::new(len, len, Self::MAX_SIZE);
        let cigar = Cigar::new(len, len);

        Self {
            read_padded,
            pattern_padded,
            matrix,
            block,
            cigar,
            len,
        }
    }

    fn resize_if_needed(&mut self, len: usize) {
        if len > self.len {
            self.read_padded = PaddedBytes::new::<NucMatrix>(len, Self::MAX_SIZE);
            self.pattern_padded = PaddedBytes::new::<NucMatrix>(len, Self::MAX_SIZE);
            self.block = Block::<true, LOCAL, LOCAL, false>::new(len, len, Self::MAX_SIZE);
            self.cigar = Cigar::new(len, len);
            self.len = len;
        }
    }
}

unsafe impl<const LOCAL: bool> Send for GlobalLocalAligner<LOCAL> {}

impl<const LOCAL: bool> Aligner for GlobalLocalAligner<LOCAL> {
    fn align(
        &mut self,
        read: &[u8],
        pattern: &[u8],
        identity_threshold: f64,
        overlap_threshold: f64,
    ) -> Option<(usize, usize, usize)> {
        self.resize_if_needed(pattern.len().max(read.len()));

        let max_size = pattern
            .len()
            .min(read.len())
            .next_power_of_two()
            .min(Self::MAX_SIZE);

        self.read_padded.set_bytes::<NucMatrix>(read, max_size);
        self.pattern_padded
            .set_bytes::<NucMatrix>(pattern, max_size);

        let min_size = if LOCAL { max_size } else { Self::MIN_SIZE };

        self.block.align(
            &self.pattern_padded,
            &self.read_padded,
            &self.matrix,
            Self::GAPS,
            min_size..=max_size,
            pattern.len() as i32,
        );

        let res = self.block.res();
        self.block.trace().cigar_eq(
            &self.pattern_padded,
            &self.read_padded,
            res.query_idx,
            res.reference_idx,
            &mut self.cigar,
        );

        let mut matches = 0;
        let mut total = 0;

        self.cigar.reverse();
        let mut read_start_idx = res.reference_idx;

        for i in 0..self.cigar.len() {
            let OpLen { op, len } = self.cigar.get(i);

            match op {
                Operation::Eq => {
                    read_start_idx -= len;
                    matches += len;
                }
                Operation::X => {
                    read_start_idx -= len;
                }
                Operation::D => {
                    read_start_idx -= len;
                }
                _ => (),
            }

            total += len;
        }

        let identity = (matches as f64) / (total as f64);
        let overlap = (matches as f64) / (pattern.len() as f64);

        if identity >= identity_threshold && overlap >= overlap_threshold {
            Some((matches, read_start_idx, res.reference_idx))
        } else {
            None
        }
    }
}

struct PrefixSuffixAligner<const PREFIX: bool> {
    read_padded: PaddedBytes,
    pattern_padded: PaddedBytes,
    matrix: NucMatrix,
    // always store trace
    block1: Block<true, true, false, true>,  // X-drop
    block2: Block<true, false, false, true>, // no X-drop
    cigar: Cigar,
    len: usize,
}

impl<const PREFIX: bool> PrefixSuffixAligner<PREFIX> {
    const MAX_SIZE: usize = 512;
    const GAPS: Gaps = Gaps {
        open: -2,
        extend: -1,
    };

    pub fn new(len: usize) -> Self {
        let read_padded = PaddedBytes::new::<NucMatrix>(len, Self::MAX_SIZE);
        let pattern_padded = PaddedBytes::new::<NucMatrix>(len, Self::MAX_SIZE);
        let matrix = NucMatrix::new_simple(1, -1);

        let block1 = Block::<true, true, false, true>::new(len, len, Self::MAX_SIZE);
        let block2 = Block::<true, false, false, true>::new(len, len, Self::MAX_SIZE);
        let cigar = Cigar::new(len, len);

        Self {
            read_padded,
            pattern_padded,
            matrix,
            block1,
            block2,
            cigar,
            len,
        }
    }

    fn resize_if_needed(&mut self, len: usize) {
        if len > self.len {
            self.read_padded = PaddedBytes::new::<NucMatrix>(len, Self::MAX_SIZE);
            self.pattern_padded = PaddedBytes::new::<NucMatrix>(len, Self::MAX_SIZE);
            self.block1 = Block::<true, true, false, true>::new(len, len, Self::MAX_SIZE);
            self.block2 = Block::<true, false, false, true>::new(len, len, Self::MAX_SIZE);
            self.cigar = Cigar::new(len, len);
            self.len = len;
        }
    }
}

unsafe impl<const PREFIX: bool> Send for PrefixSuffixAligner<PREFIX> {}

impl<const PREFIX: bool> Aligner for PrefixSuffixAligner<PREFIX> {
    fn align(
        &mut self,
        read: &[u8],
        pattern: &[u8],
        identity_threshold: f64,
        overlap_threshold: f64,
    ) -> Option<(usize, usize, usize)> {
        self.resize_if_needed(pattern.len().max(read.len()));

        let max_size = pattern
            .len()
            .min(read.len())
            .next_power_of_two()
            .min(Self::MAX_SIZE);

        if PREFIX {
            // reverse sequences to convert to aligning suffix
            self.read_padded.set_bytes_rev::<NucMatrix>(read, max_size);
            self.pattern_padded
                .set_bytes_rev::<NucMatrix>(pattern, max_size);
        } else {
            self.read_padded.set_bytes::<NucMatrix>(read, max_size);
            self.pattern_padded
                .set_bytes::<NucMatrix>(pattern, max_size);
        }

        // first align to get where the pattern starts in the read
        // note that the start gaps in the pattern are free and the alignment
        // can end whenever due to X-drop
        self.block1.align(
            &self.pattern_padded,
            &self.read_padded,
            &self.matrix,
            Self::GAPS,
            max_size..=max_size,
            pattern.len() as i32,
        );

        let res = self.block1.res();
        self.block1.trace().cigar_eq(
            &self.pattern_padded,
            &self.read_padded,
            res.query_idx,
            res.reference_idx,
            &mut self.cigar,
        );

        // use traceback to compute where the alignment started
        let mut read_start_idx = res.reference_idx;
        for i in 0..self.cigar.len() {
            let OpLen { op, len } = self.cigar.get(i);
            match op {
                Operation::Eq | Operation::X | Operation::D => read_start_idx -= len,
                _ => (),
            }
        }

        // skip second alignment if first alignment reaches the end of the read
        if res.reference_idx < read.len() {
            // get the overlapping prefix/suffix region
            if PREFIX {
                self.read_padded
                    .set_bytes::<NucMatrix>(&read[..read.len() - read_start_idx], max_size);
                self.pattern_padded
                    .set_bytes::<NucMatrix>(pattern, max_size);
            } else {
                self.read_padded
                    .set_bytes_rev::<NucMatrix>(&read[read_start_idx..], max_size);
                self.pattern_padded
                    .set_bytes_rev::<NucMatrix>(pattern, max_size);
            }

            // align again with read and pattern switched and reversed so that end gaps in the read
            // are free and the alignment ends at read_start_idx and spans the entire pattern
            self.block2.align(
                &self.read_padded,
                &self.pattern_padded,
                &self.matrix,
                Self::GAPS,
                max_size..=max_size,
                pattern.len() as i32,
            );

            let res = self.block2.res();
            self.block2.trace().cigar_eq(
                &self.read_padded,
                &self.pattern_padded,
                res.query_idx,
                res.reference_idx,
                &mut self.cigar,
            );
        }

        // count matches and total columns for calculating identity and overlap
        let mut matches = 0;
        let mut total = 0;

        for i in 0..self.cigar.len() {
            let OpLen { op, len } = self.cigar.get(i);
            if op == Operation::Eq {
                matches += len;
            }
            total += len;
        }

        let identity = (matches as f64) / (total as f64);
        let overlap = (matches as f64) / (pattern.len() as f64);

        if identity >= identity_threshold && overlap >= overlap_threshold {
            let start_idx = if PREFIX { 0 } else { read_start_idx };
            let end_idx = if PREFIX {
                read.len() - read_start_idx
            } else {
                read.len()
            };

            Some((matches, start_idx, end_idx))
        } else {
            None
        }
    }
}

#[cfg(test)]
mod edit_distance_tests {
    use super::*;

    // -- Bug 1: edit_search_myers estimates start position instead of computing it exactly --
    // Use 8bp patterns to avoid ambiguous partial matches with shorter patterns.
    #[test]
    fn test_edit_search_myers_start_position_substitution() {
        // Pattern "ACGTACGT" placed at position 10 with 1 sub in the middle (T->X)
        // text[10..18] = "ACGXACGT" vs pattern "ACGTACGT" = 1 substitution
        let text = b"NNNNNNNNNNACGXACGTNNNNNNNNNN";
        let pattern = b"ACGTACGT";
        let result = edit_search(text, pattern, 1);
        assert!(result.is_some(), "should find match with 1 edit");
        let (_matches, start, end) = result.unwrap();
        assert_eq!(end, 18, "end should be 18");
        // Buggy estimation: start = 18 - max(8-1, 8+1) = 18 - 9 = 9
        // Correct: start = 10
        assert_eq!(start, 10, "start should be exactly 10, not an estimate");
    }

    #[test]
    fn test_edit_search_myers_start_position_deletion() {
        // Pattern "ACGTACGT" with 1 deletion in text: "ACGACGT" at position 10
        // text[10..17] = "ACGACGT" aligns to "ACGTACGT" with 1 insertion (add T at pos 3)
        let text = b"NNNNNNNNNNACGACGTNNNNNNNNNN";
        let pattern = b"ACGTACGT";
        let result = edit_search(text, pattern, 1);
        assert!(result.is_some(), "should find match with 1 edit");
        let (_matches, start, end) = result.unwrap();
        // Match spans 7 text chars: text[10..17]
        assert_eq!(end, 17, "end should be 17");
        // Buggy estimation: start = 17 - max(8-1, 8+1) = 17 - 9 = 8
        // Correct: start = 10
        assert_eq!(start, 10, "start should be exactly 10, not an estimate");
    }

    // -- Bug 2: edit_prefix doesn't verify match starts at position 0 --
    #[test]
    fn test_edit_prefix_reports_correct_match_quality() {
        // text starts with NNN then has ACGT: "NNNACGTNNNN"
        // With max_edits=3, search window = 4+3 = 7: text_prefix = "NNNACGT"
        // edit_search finds exact ACGT at position 3-7 (0 edits, 4 matches)
        // But the correct PREFIX alignment is: delete NNN (3 edits), then ACGT = 1 match
        let text = b"NNNACGTNNNN";
        let pattern = b"ACGT";
        let result = edit_prefix(text, pattern, 3);
        assert!(
            result.is_some(),
            "there IS a valid prefix alignment within 3 edits"
        );
        let (matches, end) = result.unwrap();
        // Correct: prefix alignment deletes NNN (3 edits) -> matches = 4 - 3 = 1
        // Buggy: finds internal exact match (0 edits) -> matches = 4
        assert_eq!(
            matches, 1,
            "prefix match should report 1 match (3 edits for deleting NNN)"
        );
        assert_eq!(end, 7, "should consume 7 text bytes");
    }

    #[test]
    fn test_edit_prefix_with_insertion_at_start() {
        // text = "XACGT..." - 1 insertion (X) before the real prefix match
        // Correct prefix alignment: delete X (1 edit), then ACGT matches -> end=5, 1 edit
        let text = b"XACGTNNNN";
        let pattern = b"ACGT";
        let result = edit_prefix(text, pattern, 1);
        assert!(
            result.is_some(),
            "should find prefix match with 1 insertion"
        );
        let (matches, end) = result.unwrap();
        assert_eq!(matches, 3, "should report 3 matches (1 edit)");
        assert_eq!(end, 5, "should consume 5 text bytes");
    }

    // -- Bug 3: edit_search_dp estimates start position (patterns > 64bp) --
    #[test]
    fn test_edit_search_dp_start_position() {
        // 68bp pattern to force DP path (> 64bp)
        let pattern = b"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT";
        assert!(pattern.len() > 64, "pattern must be > 64bp to use DP path");
        // Place pattern at position 10 with 1 substitution in the middle (pos 34: T->N)
        let mut placed = pattern.to_vec();
        placed[34] = b'X'; // 1 substitution in the middle
                           // Use 'T' padding to be distinct from the substituted 'X'
        let mut text = vec![b'T'; 10];
        text.extend_from_slice(&placed);
        text.extend_from_slice(&[b'T'; 10]);
        let start_pos = 10;
        let end_pos = start_pos + pattern.len();

        let result = edit_search(&text, pattern, 1);
        assert!(result.is_some(), "should find match with 1 edit");
        let (_matches, start, end) = result.unwrap();
        assert_eq!(start, start_pos, "DP start should be exact");
        assert_eq!(end, end_pos, "DP end should be exact");
    }

    // -- Bug 5: Edit full match uses pattern_len as cut position --
    #[test]
    fn test_edit_distance_with_insertion() {
        // text = "ACGGT" (5bp), pattern = "ACGT" (4bp), 1 insertion (extra G)
        let text = b"ACGGT";
        let pattern = b"ACGT";
        let result = edit_distance(text, pattern, 1);
        assert!(result.is_some(), "should match with 1 edit");
    }

    #[test]
    fn test_edit_distance_with_deletion() {
        // text = "ACT" (3bp), pattern = "ACGT" (4bp), 1 deletion (missing G)
        let text = b"ACT";
        let pattern = b"ACGT";
        let result = edit_distance(text, pattern, 1);
        assert!(result.is_some(), "should match with 1 edit (deletion)");
    }

    // -- Correctness baselines --
    #[test]
    fn test_edit_distance_exact_match() {
        let result = edit_distance(b"ACGT", b"ACGT", 0);
        assert_eq!(result, Some(4));
    }

    #[test]
    fn test_edit_distance_one_sub() {
        let result = edit_distance(b"ACGC", b"ACGT", 1);
        assert_eq!(result, Some(3));
    }

    #[test]
    fn test_edit_distance_over_threshold() {
        let result = edit_distance(b"NNNN", b"ACGT", 1);
        assert_eq!(result, None);
    }

    #[test]
    fn test_edit_search_exact_match() {
        let text = b"NNNNNACGTNNNNNN";
        let pattern = b"ACGT";
        let result = edit_search(text, pattern, 0);
        assert!(result.is_some());
        let (matches, start, end) = result.unwrap();
        assert_eq!(matches, 4);
        assert_eq!(start, 5);
        assert_eq!(end, 9);
    }

    #[test]
    fn test_edit_prefix_exact() {
        let text = b"ACGTNNNNNN";
        let pattern = b"ACGT";
        let result = edit_prefix(text, pattern, 0);
        assert!(result.is_some());
        let (matches, end) = result.unwrap();
        assert_eq!(matches, 4);
        assert_eq!(end, 4);
    }

    #[test]
    fn test_edit_suffix_exact() {
        let text = b"NNNNNNACGT";
        let pattern = b"ACGT";
        let result = edit_suffix(text, pattern, 0);
        assert!(result.is_some());
        let (matches, start) = result.unwrap();
        assert_eq!(matches, 4);
        assert_eq!(start, 6);
    }

    #[test]
    fn test_edit_suffix_reports_correct_match_quality() {
        // text ends with ACGT then NNN: "NNNNACGTNNN"
        // With max_edits=3, search window covers "ACGTNNN" (last 7 bytes)
        // edit_search finds exact ACGT at position 0-4 of the window (0 edits, 4 matches)
        // But the correct SUFFIX alignment is: delete NNN at end (3 edits), ACGT matches
        let text = b"NNNNACGTNNN";
        let pattern = b"ACGT";
        let result = edit_suffix(text, pattern, 3);
        assert!(
            result.is_some(),
            "there IS a valid suffix alignment within 3 edits"
        );
        let (matches, start) = result.unwrap();
        // Correct: suffix alignment deletes trailing NNN (3 edits) -> matches = 4 - 3 = 1
        // Buggy: finds internal exact match (0 edits) -> matches = 4
        assert_eq!(
            matches, 1,
            "suffix match should report 1 match (3 edits for deleting NNN)"
        );
        assert_eq!(start, 4, "suffix match should start at position 4");
    }

    // -- Hamming distance tests --
    #[test]
    fn test_hamming_exact_match() {
        let result = hamming(b"ACGT", b"ACGT", 4);
        assert_eq!(result, Some(4));
    }

    #[test]
    fn test_hamming_one_mismatch() {
        let result = hamming(b"ACGT", b"ACGC", 3);
        assert_eq!(result, Some(3));
    }

    #[test]
    fn test_hamming_below_threshold() {
        let result = hamming(b"NNNN", b"ACGT", 4);
        assert_eq!(result, None);
    }

    #[test]
    fn test_hamming_different_lengths() {
        let result = hamming(b"ACG", b"ACGT", 3);
        assert_eq!(result, None);
    }

    #[test]
    fn test_hamming_long_sequence() {
        let a = b"ACGTACGTACGTACGT";
        let b_seq = b"ACGTACGTACGTACGT";
        let result = hamming(a, b_seq, 16);
        assert_eq!(result, Some(16));
    }

    #[test]
    fn test_hamming_long_with_mismatches() {
        let a = b"ACGTACGTACGTACGT";
        let mut b_seq = b"ACGTACGTACGTACGT".to_vec();
        b_seq[0] = b'N';
        b_seq[8] = b'N';
        let result = hamming(a, &b_seq, 14);
        assert_eq!(result, Some(14));
    }

    #[test]
    fn test_hamming_search_exact() {
        let text = b"NNNNNACGTNNNNNN";
        let pattern = b"ACGT";
        let result = hamming_search(text, pattern, 4);
        assert!(result.is_some());
        let (matches, start, end) = result.unwrap();
        assert_eq!(matches, 4);
        assert_eq!(start, 5);
        assert_eq!(end, 9);
    }

    #[test]
    fn test_hamming_search_with_mismatch() {
        let text = b"NNNNNACGCNNNNNN";
        let pattern = b"ACGT";
        let result = hamming_search(text, pattern, 3);
        assert!(result.is_some());
        let (matches, start, end) = result.unwrap();
        assert_eq!(matches, 3);
        assert_eq!(start, 5);
        assert_eq!(end, 9);
    }

    #[test]
    fn test_hamming_search_no_match() {
        let text = b"NNNNNNNNNNNN";
        let pattern = b"ACGT";
        let result = hamming_search(text, pattern, 4);
        assert!(result.is_none());
    }

    // -- Edit distance edge cases --
    #[test]
    fn test_edit_distance_empty_pattern() {
        let result = edit_distance(b"ACGT", b"", 0);
        assert_eq!(result, None);
    }

    #[test]
    fn test_edit_distance_empty_text() {
        let result = edit_distance(b"", b"ACGT", 4);
        assert_eq!(result, Some(0));
    }

    #[test]
    fn test_edit_distance_both_empty() {
        let result = edit_distance(b"", b"", 0);
        assert_eq!(result, Some(0));
    }

    #[test]
    fn test_edit_distance_length_diff_exceeds_max() {
        let result = edit_distance(b"A", b"ACGTACGT", 2);
        assert_eq!(result, None);
    }

    #[test]
    fn test_edit_search_empty_pattern() {
        let result = edit_search(b"ACGT", b"", 0);
        // Empty pattern behavior depends on implementation
        // Just verify it does not panic
        let _ = result;
    }

    #[test]
    fn test_edit_prefix_empty_pattern() {
        let result = edit_prefix(b"ACGT", b"", 0);
        assert!(result.is_some());
        assert_eq!(result.unwrap(), (0, 0));
    }

    #[test]
    fn test_edit_suffix_empty_pattern() {
        let result = edit_suffix(b"ACGT", b"", 0);
        assert!(result.is_some());
    }

    #[test]
    fn test_edit_suffix_with_insertion_at_end() {
        // text = "NNNNACGTX" - 1 insertion (X) after the real suffix match
        // Correct suffix alignment: delete X (1 edit), then ACGT matches -> start=4, 1 edit
        let text = b"NNNNACGTX";
        let pattern = b"ACGT";
        let result = edit_suffix(text, pattern, 1);
        assert!(
            result.is_some(),
            "should find suffix match with 1 insertion"
        );
        let (matches, start) = result.unwrap();
        assert_eq!(matches, 3, "should report 3 matches (1 edit)");
        assert_eq!(start, 4, "should start at position 4");
    }
}
