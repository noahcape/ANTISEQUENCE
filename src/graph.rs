use std::marker::{Send, Sync};
use std::ops::RangeBounds;
use std::path::Path;
use std::sync::Arc;
use std::thread;
use std::sync::OnceLock;

use crate::errors::*;
use crate::expr::*;
use crate::read::*;
use crate::trace::*;

mod ops;
pub use ops::*;

/// Computation graph of read operations, where each operation is a node.
pub struct Graph<T: Trace = NoTrace> {
    nodes: Vec<Arc<dyn GraphNode<T>>>,
}

#[derive(Debug, Clone)]
pub struct MatchDistanceCounts {
    pub label: String,
    pub counts: Vec<usize>,
    pub total: usize,
}

#[derive(Debug, Clone)]
pub struct InputStats {
    pub n_fastqs: usize,
    pub read_counts: Vec<usize>,
    pub read_length_min: Vec<usize>,
    pub read_length_max: Vec<usize>,
    pub read_length_sum: Vec<usize>,
}

pub trait GraphNode<T: Trace = NoTrace>: Send + Sync {
    #[inline(always)]
    fn run(&self, reads: Option<Vec<Read>>, trace: &T) -> Result<(Option<Vec<Read>>, bool)> {
        let start = trace.start(&reads);
        let Some(reads) = reads else {
            panic!("Expected some reads!")
        };
        let res = self.run_inner(reads)?;
        trace.add(self.name(), start, &res.0);
        Ok(res)
    }
    fn run_inner(&self, _reads: Vec<Read>) -> Result<(Option<Vec<Read>>, bool)> {
        unimplemented!()
    }
    fn required_names(&self) -> &[LabelOrAttr];
    fn name(&self) -> &'static str;

    /// Optional hook for nodes that expose match distance statistics.
    ///
    /// Default implementation returns `None` so that most nodes do not need
    /// to be aware of statistics collection.
    #[inline]
    fn match_distance_counts(&self) -> Option<MatchDistanceCounts> {
        None
    }

    #[inline]
    fn input_stats(&self) -> Option<InputStats> {
        None
    }
}

impl<T: Trace> Graph<T> {
    /// Create a new empty graph.
    pub fn new() -> Self {
        Self { nodes: Vec::new() }
    }

    /// Add a read operation node to the graph and return the node.
    pub fn add<G: GraphNode<T> + 'static>(&mut self, node: G) -> Arc<G> {
        let a = Arc::new(node);
        let b = Arc::clone(&a);
        self.nodes.push(a);
        b
    }

    /// Run a graph until all reads processed.
    pub fn run(&self) -> Result<()> {
        self.run_trace(DEFAULT_TRACE_PATH)
    }

    /// Run a graph until all reads processed, outputting the trace to the specified path.
    pub fn run_trace(&self, trace_path: impl AsRef<Path>) -> Result<()> {
        let trace = T::new(trace_path);
        let res = self.run_trace_inner(&trace);
        trace.finish();
        res
    }

    fn run_trace_inner(&self, trace: &T) -> Result<()> {
        let mut next_input: Option<Vec<Read>> = None;
        loop {
            // Pass next_input to recycle the vector
            let (out, done) = self.run_one(next_input, trace)?;

            // Recycle output vector for next input, but DO NOT clear.
            // We let InputFastqOp handle the clearing/recycling logic to reuse Read internal buffers.
            next_input = out;

            if done {
                break;
            }
        }

        Ok(())
    }

    /// Run a graph in parallel (multithreading) until all reads processed.
    pub fn run_with_threads(&self, threads: usize) {
        self.run_with_threads_trace(threads, DEFAULT_TRACE_PATH);
    }

    /// Run a graph in parallel (multithreading) until all reads processed, with tracing.
    pub fn run_with_threads_trace(&self, threads: usize, trace_path: impl AsRef<Path>) {
        let trace = T::new(trace_path);
        self.run_with_threads_trace_inner(threads, &trace);
        trace.finish();
    }

    fn run_with_threads_trace_inner(&self, threads: usize, trace: &T) {
        assert!(threads >= 1, "Number of threads must be greater than zero");

        thread::scope(|s| {
            for _ in 0..threads {
                s.spawn(move || {
                    self.run_trace_inner(&trace)
                        .unwrap_or_else(|e| panic!("{e}"));
                });
            }
        });
    }

    /// Collect per-node match distance histograms from all nodes that expose
    /// them via `GraphNode::match_distance_counts`.
    ///
    /// Each entry corresponds to a single node instance and contains the label
    /// name (e.g. "seq1.brc") and a vector of counts indexed by edit
    /// distance (0, 1, 2, ...).
    pub fn match_distance_counts(&self) -> Vec<MatchDistanceCounts> {
        let mut out = Vec::new();
        for node in &self.nodes {
            if let Some(counts) = node.match_distance_counts() {
                out.push(counts);
            }
        }
        out
    }

    /// Collect input statistics from the first node that exposes them
    /// via `GraphNode::input_stats` (typically the InputFastqOp).
    pub fn input_stats(&self) -> Option<InputStats> {
        for node in &self.nodes {
            if let Some(stats) = node.input_stats() {
                return Some(stats);
            }
        }
        None
    }

    /// Run a single batch of reads through the graph.
    ///
    /// Returns an additional boolean indicating whether the graph is done executing.
    /// If the required label or attribute names for an operation are not available,
    /// the the operation is skipped.
    #[inline(always)]
    pub fn run_one(&self, mut curr: Option<Vec<Read>>, trace: &T) -> Result<(Option<Vec<Read>>, bool)> {
        let trust = trust_required_checks();
        for node in &self.nodes {
            // If there is no current read, only the input node can produce one.
            if curr.is_none() {
                let (c, done) = node.run(None, trace)?;
                curr = c;
                if done { return Ok((curr, done)); }
                if curr.is_none() { break; }
                continue;
            }

            // Skip nodes whose requirements are not satisfied, unless trusted.
            // Heuristic: Check the first read as a representative.
            if !trust && !node.required_names().is_empty() {
                if let Some(reads) = &curr {
                    if let Some(first) = reads.first() {
                        if !first.has_names(node.required_names()) {
                            continue;
                        }
                    }
                }
            }

            // Call node.run so nodes that override run (and not run_inner) still work.
            let (c, done) = node.run(curr, trace)?;
            curr = c;

            if done { return Ok((curr, done)); }
            if curr.is_none() { break; }
        }

        Ok((curr, false))
    }

    /// Try running a single batch of reads through the graph.
    ///
    /// Returns two booleans: the first one is whether the read has "failed" (does not have
    /// a required label or attribute name) and the second one is whether the graph is done
    /// executing.
    pub fn try_run_one(
        &self,
        mut curr: Option<Vec<Read>>,
        trace: &T,
    ) -> Result<(Option<Vec<Read>>, bool, bool)> {
        for node in &self.nodes {
            if let Some(reads) = &curr {
                if let Some(first) = reads.first() {
                     if !first.has_names(node.required_names()) {
                        return Ok((curr, true, false));
                    }
                }
            }

            let (c, done) = node.run(curr, trace)?;
            curr = c;

            if done {
                return Ok((curr, false, done));
            }
            if curr.is_none() {
                break;
            }
        }

        Ok((curr, false, false))
    }
}

#[inline(always)]
fn trust_required_checks() -> bool {
    static TRUST: OnceLock<bool> = OnceLock::new();
    *TRUST.get_or_init(|| {
        std::env::var("ANTISEQ_TRUST_NAMES")
            .ok()
            .map(|v| v == "1" || v.eq_ignore_ascii_case("true"))
            .unwrap_or(false)
    })
}

pub use MatchType::*;
pub use Threshold::*;

/// Algorithm types for matching patterns.
///
/// For alignment-based algorithms, `sequence identity = matches / (matches + mismatches + insertions + deletions)`
/// and `overlap = matches / pattern_length`.
///
/// Insertions and deletions that are not part of the alignment are not included in the sequence
/// identity computation. This is important for local alignment, where the start and end of the
/// pattern can be excluded from the alignment, and prefix/suffix alignment, where the start/end
/// of the pattern can be excluded from the alignment (prefix/suffix "overhang").
#[derive(Copy, Clone, Debug, PartialEq)]
pub enum MatchType {
    /// Exact match.
    ///
    /// A match will result in one new interval: the entire string.
    Exact,
    /// Exact prefix match.
    ///
    /// A match will result in two new interval: the matched prefix and the rest of the string.
    ExactPrefix,
    /// Exact suffix match.
    ///
    /// A match will result in two new interval: the rest of the string and the matched
    /// suffix.
    ExactSuffix,
    /// Exact match search.
    ///
    /// A match will result in three new interval: everything before the exact match, the exact
    /// matching region, everything after the exact match.
    ExactSearch,
    /// Hamming-distance-based matching.
    ///
    /// Threshold is for the number of matching bases.
    ///
    /// A match will result in one new interval: the entire string.
    Hamming(Threshold),
    /// Hamming-distance-based prefix matching.
    ///
    /// Threshold is for the number of matching bases.
    ///
    /// A match will result in two new interval: the matched prefix and the rest of the
    /// string.
    HammingPrefix(Threshold),
    /// Hamming-distance-based suffix matching.
    ///
    /// Threshold is for the number of matching bases.
    ///
    /// A match will result in two new interval: the rest of the string and the matched
    /// suffix.
    HammingSuffix(Threshold),
    /// Hamming-distance-based searching.
    ///
    /// Threshold is for the number of matching bases.
    ///
    /// A match will result in three new interval: everything before the match, the matching
    /// region, and everything after the match.
    HammingSearch(Threshold),
    /// Global-alignment-based matching.
    ///
    /// Threshold is for the sequence identity.
    ///
    /// A match will result in one new interval: the entire string.
    GlobalAln(f64),
    /// Local-alignment-based matching.
    ///
    /// A match will result in three new interval: everything before the aligned region, the locally aligned
    /// region, and everything after the aligned region.
    LocalAln { identity: f64, overlap: f64 },
    /// Prefix-alignment-based matching.
    ///
    /// A match will result in two new interval: the matched prefix and the rest of the
    /// string.
    PrefixAln { identity: f64, overlap: f64 },
    /// Suffix-alignment-based matching.
    ///
    /// A match will result in two new interval: the rest of the string and the matched
    /// suffix.
    SuffixAln { identity: f64, overlap: f64 },
    /// Exact-alignment within a range.
    ///
    /// A match will result in three new intervals: everything before the aligned region, the aligned
    /// region, and everything after the aligned region.
    /// Use inclusive range indexing, from..=to
    ExactBoundedMatch { from: usize, to: usize },
    /// Hamming-distance-based alignment within a range.
    ///
    /// A match will result in three new intervals: everything before the aligned region, the aligned
    /// region, and everything after the aligned region.
    /// Use inclusive range indexing, from..=to
    HammingBoundedMatch {
        threshold: Threshold,
        from: usize,
        to: usize,
    },
}

impl MatchType {
    pub fn num_mappings(&self) -> usize {
        use MatchType::*;
        match self {
            Exact | Hamming(_) | GlobalAln(_) => 1,
            ExactPrefix
            | ExactSuffix
            | HammingPrefix(_)
            | HammingSuffix(_)
            | PrefixAln { .. }
            | SuffixAln { .. } => 2,
            ExactSearch
            | HammingSearch(_)
            | LocalAln { .. }
            | HammingBoundedMatch { .. }
            | ExactBoundedMatch { .. } => 3,
        }
    }

    pub fn k(&self, len: usize) -> usize {
        let k_from_edits = |len: usize, e: usize| (len - e).div_ceil(e + 1);
        use MatchType::*;
        match self {
            Exact => len,
            ExactPrefix => len,
            ExactSuffix => len,
            ExactSearch => len,
            ExactBoundedMatch { .. } => len,
            Hamming(t) => k_from_edits(len, t.get(len)),
            HammingPrefix(t) => k_from_edits(len, t.get(len)),
            HammingSuffix(t) => k_from_edits(len, t.get(len)),
            HammingSearch(t) => k_from_edits(len, t.get(len)),
            HammingBoundedMatch { threshold: t, .. } => k_from_edits(len, t.get(len)),
            GlobalAln(identity) => {
                k_from_edits(len, len - (((len as f64) * identity).ceil() as usize))
            }
            PrefixAln { identity, overlap } => {
                let len = ((len as f64) * overlap).ceil() as usize;
                k_from_edits(len, len - (((len as f64) * identity).ceil() as usize))
            }
            SuffixAln { identity, overlap } => {
                let len = ((len as f64) * overlap).ceil() as usize;
                k_from_edits(len, len - (((len as f64) * identity).ceil() as usize))
            }
            LocalAln { identity, overlap } => {
                let len = ((len as f64) * overlap).ceil() as usize;
                k_from_edits(len, len - (((len as f64) * identity).ceil() as usize))
            }
        }
    }
}

/// Either a count or a fraction.
///
/// Typically used for specifying the similarity threshold when matching patterns.
/// The fraction is typically of the length of the pattern.
#[derive(Copy, Clone, Debug, PartialEq)]
pub enum Threshold {
    Count(usize),
    Frac(f64),
}

impl Threshold {
    pub fn get(&self, len: usize) -> usize {
        use Threshold::*;
        match self {
            Count(c) => *c,
            Frac(f) => (*f * (len as f64)) as usize,
        }
    }
}
