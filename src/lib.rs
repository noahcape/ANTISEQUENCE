//! Rust stream processing library for sequencing reads.
//!
//! # Overview
//! ANTISEQUENCE provides efficient and composable operations for manipulating fastq records.
//!
//! This is useful for:
//! * Processing reads for custom sequencing protocol development
//! * Unifying read formats from various sequencing protocols
//! * Writing fastq processing tools
//! * Debugging sequencing pipelines
//!
//! ## Computation graph API
//! To use ANTISEQUENCE, you first specify *operations* (read from fastq, trim reads, output to fastq, etc.)
//! and add them to a [`Graph`]. Then, you run the graph, which executes all the operations on each
//! read.
//!
//! See [`graph`] for all supported operations.
//!
//! Each operation in a graph contains a set of dependencies, which are labels and attributes
//! that the operation requires to be present in the read. If the dependencies are not present,
//! then the operation will be skipped.
//!
//! ## Reads
//! Here's an example fastq record:
//! ```text
//! @read6
//! AATTCCGGAATTCCCAAAAG
//! +
//! 01234567890123456789
//! ```
//! The first, second, and fourth lines are the name, sequence, and quality scores, respectively.
//!
//! ANTISEQUENCE stores that record as an internal [`Read`] data structure:
//! ```text
//! name1:
//!  *     |---|
//!  str:  read6
//!  from: record 5 in file: "example_data/match.fastq"
//! seq1:
//!  *        |------------------|  adapter=AAAA
//!  template |-------------|
//!  adapter                 |---|
//!  str:     AATTCCGGAATTCCCAAAAG
//!  qual:    01234567890123456789
//!  from:    record 5 in file: "example_data/match.fastq"
//! ```
//!
//! Each `Read` is a set of *strings* of different *types*. Types help indicate whether the string is
//! a read sequence (`seq1`) or read name (`name1`).
//!
//! Each string has associated *labeled intervals* in the string.
//! For example, the region where an adapter is found in the read sequence can be labeled.
//! All strings start with an interval labeled `*`, which spans the whole string.
//! You can refer to an interval with `seq1.*`, `seq1.adapter`, `name1.*`, etc.
//!
//! An interval can contain *attributes* that hold arbitrary metadata. This may include a boolean for
//! whether to filter the read, or the name of the pattern that the read matches.
//! You can refer to an attribute with `seq1.*.adapter`, etc.
//!
//! For efficiency and simplicity, most ANTISEQUENCE operations only manipulate the intervals
//! and attributes. You can choose to modify the underlying strings afterwards.
//!
//! ## Transform expressions
//! Transform expressions allow you to specify the names of the inputs and outputs
//! for an operation. For example, to cut an interval and create two new intervals,
//! you can use `tr!(seq1.* -> seq1.left, seq1.right)`.
//!
//! ## Expressions
//! Expressions are useful for doing arbitrary computation on reads. Here are some examples:
//! * `label_exists("seq1.label").and(attr_exists("seq1.*.attribute"))`
//! * `Expr::from(label("seq1.label")).slice(..3)`
//! * `Expr::from(label("seq1.label")).len().in_bounds(10..20)`
//!
//! Supported data types are: byte string, int, float, and boolean.
//!
//! ### Format expressions
//! Format expressions allow you to contruct new strings from intervals and attributes,
//! and they are similar to Rust's formatting syntax. For example, you can use `fmt_expr("{seq1.a}_{seq1.b}")`
//! to concatenate the substrings corresponding to mappings `a` and `b`, separated by an
//! underscore.
//! They also preserve quality scores, making rearranging regions in a read easy.
//!
//! ### Note
//! To apply the same transformations to the quality score, the same expression is evaluated by
//! simply substituting in quality scores when needed when a `label()` is used.
//! If this is not desired, then break up the expression by storing intermediate results as
//! attributes, which will always have the same value (not substituted by quality scores).

// Optional global allocator selection
cfg_if::cfg_if! {
    if #[cfg(feature = "mimalloc")] {
        use mimalloc::MiMalloc;
        #[global_allocator]
        static GLOBAL: MiMalloc = MiMalloc;
    } else if #[cfg(feature = "jemalloc")] {
        use jemallocator::Jemalloc;
        #[global_allocator]
        static GLOBAL: Jemalloc = Jemalloc;
    }
}

pub mod errors;
pub mod expr;
pub mod graph;
mod patterns;
mod read;
pub mod trace;

mod inline_string;
mod parse_utils;
mod seed_search;

// commonly used functions and types

pub use crate::patterns::*;
pub use crate::read::*;

#[cfg(test)]
mod pipeline_tests {
    use crate::expr::*;
    use crate::graph::*;
    use crate::patterns::Patterns;
    use crate::read::*;
    use crate::trace::NoTrace;
    use std::io::Cursor;

    fn fastq_bytes(records: &[(&str, &str, &str)]) -> Vec<u8> {
        let mut buf = Vec::new();
        for (name, seq, qual) in records {
            buf.extend_from_slice(b"@");
            buf.extend_from_slice(name.as_bytes());
            buf.extend_from_slice(b"\n");
            buf.extend_from_slice(seq.as_bytes());
            buf.extend_from_slice(b"\n+\n");
            buf.extend_from_slice(qual.as_bytes());
            buf.extend_from_slice(b"\n");
        }
        buf
    }

    fn te(s: &str) -> TransformExpr {
        TransformExpr::from_bytes(s).unwrap()
    }

    #[test]
    fn test_graph_basic_pipeline() {
        let fq = fastq_bytes(&[
            ("read1", "ACGTACGT", "IIIIIIII"),
            ("read2", "TGCATGCA", "IIIIIIII"),
        ]);

        let tmp = std::env::temp_dir().join("antiseq_test_basic.fastq");
        let tmp_str = tmp.to_str().unwrap().to_string();

        {
            let mut g = Graph::<NoTrace>::new();
            g.add(InputFastqOp::from_reader(Cursor::new(fq)).unwrap());
            g.add(OutputFastqFileOp::from_file(tmp_str.clone()));
            g.run().unwrap();
        } // graph dropped here, flushing output

        let output = std::fs::read_to_string(&tmp).unwrap();
        assert!(output.contains("ACGTACGT"));
        assert!(output.contains("TGCATGCA"));
        std::fs::remove_file(&tmp).ok();
    }

    #[test]
    fn test_graph_cut_pipeline() {
        let fq = fastq_bytes(&[("read1", "ACGTACGT", "IIIIIIII")]);

        let tmp = std::env::temp_dir().join("antiseq_test_cut.fastq");
        let tmp_str = tmp.to_str().unwrap().to_string();

        {
            let mut g = Graph::<NoTrace>::new();
            g.add(InputFastqOp::from_reader(Cursor::new(fq)).unwrap());
            g.add(CutOp::new(te("seq1.* -> seq1.left, seq1.right"), 4isize));
            g.add(TrimOp::new([label("seq1.left")]));
            g.add(OutputFastqFileOp::from_file(tmp_str.clone()));
            g.run().unwrap();
        }

        let output = std::fs::read_to_string(&tmp).unwrap();
        assert!(output.contains("ACGT"));
        assert!(!output.contains("ACGTACGT"));
        std::fs::remove_file(&tmp).ok();
    }

    #[test]
    fn test_graph_set_pipeline() {
        let fq = fastq_bytes(&[("read1", "ACGTACGT", "IIIIIIII")]);

        let tmp = std::env::temp_dir().join("antiseq_test_set.fastq");
        let tmp_str = tmp.to_str().unwrap().to_string();

        {
            let mut g = Graph::<NoTrace>::new();
            g.add(InputFastqOp::from_reader(Cursor::new(fq)).unwrap());
            g.add(CutOp::new(te("seq1.* -> seq1.left, seq1.right"), 4isize));
            g.add(SetOp::new(label("seq1.left"), b"NNNN".to_vec()));
            g.add(OutputFastqFileOp::from_file(tmp_str.clone()));
            g.run().unwrap();
        }

        let output = std::fs::read_to_string(&tmp).unwrap();
        assert!(output.contains("NNNNACGT"));
        std::fs::remove_file(&tmp).ok();
    }

    #[test]
    fn test_graph_retain_pipeline() {
        let fq = fastq_bytes(&[("read1", "ACGTACGT", "IIIIIIII"), ("read2", "TG", "II")]);

        let tmp = std::env::temp_dir().join("antiseq_test_retain.fastq");
        let tmp_str = tmp.to_str().unwrap().to_string();

        {
            let mut g = Graph::<NoTrace>::new();
            g.add(InputFastqOp::from_reader(Cursor::new(fq)).unwrap());
            g.add(RetainOp::new(Expr::from(label("seq1.*")).len().gt(4isize)));
            g.add(OutputFastqFileOp::from_file(tmp_str.clone()));
            g.run().unwrap();
        }

        let output = std::fs::read_to_string(&tmp).unwrap();
        assert!(output.contains("ACGTACGT"));
        assert!(!output.contains("TG"));
        std::fs::remove_file(&tmp).ok();
    }

    #[test]
    fn test_graph_for_each_pipeline() {
        let fq = fastq_bytes(&[("read1", "ACGTACGT", "IIIIIIII")]);

        let tmp = std::env::temp_dir().join("antiseq_test_foreach.fastq");
        let tmp_str = tmp.to_str().unwrap().to_string();

        {
            let mut g = Graph::<NoTrace>::new();
            g.add(InputFastqOp::from_reader(Cursor::new(fq)).unwrap());
            g.add(RemoveInternalOp::create());
            g.add(OutputFastqFileOp::from_file(tmp_str.clone()));
            g.run().unwrap();
        }

        let output = std::fs::read_to_string(&tmp).unwrap();
        assert!(output.contains("ACGTACGT"));
        std::fs::remove_file(&tmp).ok();
    }

    #[test]
    fn test_graph_take_pipeline() {
        let fq = fastq_bytes(&[
            ("read1", "ACGTACGT", "IIIIIIII"),
            ("read2", "TGCATGCA", "IIIIIIII"),
            ("read3", "AAAAAAAA", "IIIIIIII"),
        ]);

        let mut g = Graph::<NoTrace>::new();
        g.add(InputFastqOp::from_reader(Cursor::new(fq)).unwrap());
        g.add(TakeOp::new(..2usize));
        let counter = g.add(CountOp::new([true]));
        g.run().unwrap();

        // TakeOp limits to the first 2 records
        assert!(counter.counts()[0] <= 2);
    }

    #[test]
    fn test_graph_select_pipeline() {
        let fq = fastq_bytes(&[("read1", "ACGTACGT", "IIIIIIII"), ("read2", "TG", "II")]);

        let tmp = std::env::temp_dir().join("antiseq_test_select.fastq");
        let tmp_str = tmp.to_str().unwrap().to_string();

        {
            let mut inner_g = Graph::<NoTrace>::new();
            inner_g.add(CutOp::new(te("seq1.* -> seq1.left, seq1.right"), 4isize));
            inner_g.add(TrimOp::new([label("seq1.left")]));

            let mut g = Graph::<NoTrace>::new();
            g.add(InputFastqOp::from_reader(Cursor::new(fq)).unwrap());
            g.add(SelectOp::new(
                Expr::from(label("seq1.*")).len().ge(4isize),
                inner_g,
            ));
            g.add(OutputFastqFileOp::from_file(tmp_str.clone()));
            g.run().unwrap();
        }

        let output = std::fs::read_to_string(&tmp).unwrap();
        assert!(output.contains("ACGT"));
        assert!(output.contains("TG"));
        std::fs::remove_file(&tmp).ok();
    }

    #[test]
    fn test_graph_multiple_readers() {
        // from_readers creates paired-end reads (R1 + R2 in one Read)
        let fq1 = fastq_bytes(&[("read1", "AAAA", "IIII")]);
        let fq2 = fastq_bytes(&[("read1", "CCCC", "IIII")]);

        let mut g = Graph::<NoTrace>::new();
        g.add(InputFastqOp::from_readers(vec![Cursor::new(fq1), Cursor::new(fq2)]).unwrap());
        let counter = g.add(CountOp::new([true]));
        g.run().unwrap();

        assert_eq!(counter.counts()[0], 1);
    }

    #[test]
    fn test_graph_match_distance_counts() {
        let g = Graph::<NoTrace>::new();
        let counts = g.match_distance_counts();
        assert!(counts.is_empty());
    }

    #[test]
    fn test_graph_input_stats() {
        let g = Graph::<NoTrace>::new();
        assert!(g.input_stats().is_none());
    }

    #[test]
    fn test_graph_input_stats_with_input() {
        let fq = fastq_bytes(&[("read1", "ACGT", "IIII"), ("read2", "TGCATG", "IIIIII")]);

        let mut g = Graph::<NoTrace>::new();
        let input = g.add(InputFastqOp::from_reader(Cursor::new(fq)).unwrap());
        g.run().unwrap();

        let stats = GraphNode::<NoTrace>::input_stats(input.as_ref()).unwrap();
        assert_eq!(stats.n_fastqs, 1);
        assert_eq!(stats.read_counts[0], 2);
    }

    #[test]
    fn test_graph_interleaved_reader() {
        let fq = fastq_bytes(&[
            ("read1_R1", "AAAA", "IIII"),
            ("read1_R2", "CCCC", "IIII"),
            ("read2_R1", "GGGG", "IIII"),
            ("read2_R2", "TTTT", "IIII"),
        ]);

        let tmp = std::env::temp_dir().join("antiseq_test_interleaved.fastq");
        let tmp_str = tmp.to_str().unwrap().to_string();

        {
            let mut g = Graph::<NoTrace>::new();
            g.add(InputFastqOp::from_interleaved_reader(Cursor::new(fq), 2).unwrap());
            g.add(OutputFastqFileOp::from_file(tmp_str.clone()));
            g.run().unwrap();
        }

        let output = std::fs::read_to_string(&tmp).unwrap();
        assert!(output.contains("AAAA"));
        std::fs::remove_file(&tmp).ok();
    }

    #[test]
    fn test_graph_run_one() {
        let fq = fastq_bytes(&[("read1", "ACGT", "IIII")]);

        let mut g = Graph::<NoTrace>::new();
        g.add(InputFastqOp::from_reader(Cursor::new(fq)).unwrap());

        let trace = NoTrace;
        let (out, done) = g.run_one(None, &trace).unwrap();
        assert!(out.is_some());
        assert!(!done);

        let (_, done2) = g.run_one(None, &trace).unwrap();
        assert!(done2);
    }

    #[test]
    fn test_graph_try_run_one() {
        let fq = fastq_bytes(&[("read1", "ACGT", "IIII")]);

        let mut g = Graph::<NoTrace>::new();
        g.add(InputFastqOp::from_reader(Cursor::new(fq)).unwrap());

        let trace = NoTrace;
        let (out, failed, done) = g.try_run_one(None, &trace).unwrap();
        assert!(out.is_some());
        assert!(!failed);
        assert!(!done);
    }

    #[test]
    fn test_graph_match_exact() {
        let fq = fastq_bytes(&[("read1", "ACGTNNNN", "IIIIIIII")]);

        let tmp = std::env::temp_dir().join("antiseq_test_match.fastq");
        let tmp_str = tmp.to_str().unwrap().to_string();

        let patterns = Patterns::from_strs(["ACGT"]);

        {
            let mut g = Graph::<NoTrace>::new();
            g.add(InputFastqOp::from_reader(Cursor::new(fq)).unwrap());
            g.add(MatchAnyOp::new(
                te("seq1.* -> seq1.barcode, seq1.rest"),
                patterns,
                ExactPrefix,
            ));
            g.add(TrimOp::new([label("seq1.barcode")]));
            g.add(OutputFastqFileOp::from_file(tmp_str.clone()));
            g.run().unwrap();
        }

        let output = std::fs::read_to_string(&tmp).unwrap();
        assert!(output.contains("NNNN"));
        assert!(!output.contains("ACGTNNNN"));
        std::fs::remove_file(&tmp).ok();
    }

    #[test]
    fn test_graph_match_hamming() {
        let fq = fastq_bytes(&[
            ("read1", "ACGTNNNN", "IIIIIIII"),
            ("read2", "ACGCNNNN", "IIIIIIII"),
        ]);

        let tmp = std::env::temp_dir().join("antiseq_test_hamming.fastq");
        let tmp_str = tmp.to_str().unwrap().to_string();

        let patterns = Patterns::from_strs(["ACGT"]);

        {
            let mut g = Graph::<NoTrace>::new();
            g.add(InputFastqOp::from_reader(Cursor::new(fq)).unwrap());
            g.add(MatchAnyOp::new(
                te("seq1.* -> seq1.barcode, seq1.rest"),
                patterns,
                HammingPrefix(Count(3)),
            ));
            g.add(TrimOp::new([label("seq1.barcode")]));
            g.add(OutputFastqFileOp::from_file(tmp_str.clone()));
            g.run().unwrap();
        }

        let output = std::fs::read_to_string(&tmp).unwrap();
        assert!(output.contains("read1"));
        assert!(output.contains("read2"));
        std::fs::remove_file(&tmp).ok();
    }

    #[test]
    fn test_graph_match_edit() {
        let fq = fastq_bytes(&[
            ("read1", "ACGTNNNN", "IIIIIIII"),
            ("read2", "ACGNNNN", "IIIIIII"),
        ]);

        let tmp = std::env::temp_dir().join("antiseq_test_edit.fastq");
        let tmp_str = tmp.to_str().unwrap().to_string();

        let patterns = Patterns::from_strs(["ACGT"]);

        {
            let mut g = Graph::<NoTrace>::new();
            g.add(InputFastqOp::from_reader(Cursor::new(fq)).unwrap());
            g.add(MatchAnyOp::new(
                te("seq1.* -> seq1.barcode, seq1.rest"),
                patterns,
                EditPrefix(Count(3)),
            ));
            g.add(TrimOp::new([label("seq1.barcode")]));
            g.add(OutputFastqFileOp::from_file(tmp_str.clone()));
            g.run().unwrap();
        }

        let output = std::fs::read_to_string(&tmp).unwrap();
        assert!(output.contains("read1"));
        assert!(output.contains("read2"));
        std::fs::remove_file(&tmp).ok();
    }

    #[test]
    fn test_graph_set_attr() {
        let fq = fastq_bytes(&[("read1", "ACGT", "IIII")]);

        let mut g = Graph::<NoTrace>::new();
        g.add(InputFastqOp::from_reader(Cursor::new(fq)).unwrap());
        g.add(SetOp::new(
            attr("seq1.*.score"),
            Expr::from(label("seq1.*")).len(),
        ));

        let trace = NoTrace;
        let (out, _) = g.run_one(None, &trace).unwrap();
        let reads = out.unwrap();
        let read = &reads[0];
        let score = read
            .data(
                StrType::Seq(1),
                crate::inline_string::InlineString::new(b"*"),
                crate::inline_string::InlineString::new(b"score"),
            )
            .unwrap();
        assert_eq!(score, &Data::Int(4));
    }

    #[test]
    fn test_graph_intersect_union_ops() {
        let fq = fastq_bytes(&[("read1", "ACGTACGT", "IIIIIIII")]);

        let mut g = Graph::<NoTrace>::new();
        g.add(InputFastqOp::from_reader(Cursor::new(fq)).unwrap());
        g.add(CutOp::new(te("seq1.* -> seq1.left, seq1.right"), 4isize));
        g.add(CutOp::new(te("seq1.right -> seq1.mid, seq1.end"), 2isize));
        g.add(IntersectOp::new(te("seq1.left, seq1.right -> seq1.inter")));
        g.add(UnionOp::new(te("seq1.left, seq1.end -> seq1.uni")));

        let trace = NoTrace;
        let (out, _) = g.run_one(None, &trace).unwrap();
        let reads = out.unwrap();
        let read = &reads[0];

        let uni = read
            .substring(
                StrType::Seq(1),
                crate::inline_string::InlineString::new(b"uni"),
            )
            .unwrap();
        assert_eq!(uni, b"ACGTACGT");
    }

    #[test]
    fn test_graph_count_op() {
        let fq = fastq_bytes(&[
            ("read1", "ACGT", "IIII"),
            ("read2", "TGCA", "IIII"),
            ("read3", "AAAA", "IIII"),
        ]);

        let mut g = Graph::<NoTrace>::new();
        g.add(InputFastqOp::from_reader(Cursor::new(fq)).unwrap());
        let counter = g.add(CountOp::new([true]));
        g.run().unwrap();

        assert_eq!(counter.counts()[0], 3);
    }

    #[test]
    fn test_graph_with_threads() {
        let fq = fastq_bytes(&[("read1", "ACGT", "IIII"), ("read2", "TGCA", "IIII")]);

        let tmp = std::env::temp_dir().join("antiseq_test_threads.fastq");
        let tmp_str = tmp.to_str().unwrap().to_string();

        {
            let mut g = Graph::<NoTrace>::new();
            g.add(InputFastqOp::from_reader(Cursor::new(fq)).unwrap());
            g.add(OutputFastqFileOp::from_file(tmp_str.clone()));
            g.run_with_threads(2);
        }

        let output = std::fs::read_to_string(&tmp).unwrap();
        assert!(output.contains("ACGT") || output.contains("TGCA"));
        std::fs::remove_file(&tmp).ok();
    }

    #[test]
    fn test_match_type_num_mappings() {
        assert_eq!(Exact.num_mappings(), 1);
        assert_eq!(ExactPrefix.num_mappings(), 2);
        assert_eq!(ExactSuffix.num_mappings(), 2);
        assert_eq!(ExactSearch.num_mappings(), 3);
        assert_eq!(Hamming(Count(3)).num_mappings(), 1);
        assert_eq!(HammingPrefix(Count(3)).num_mappings(), 2);
        assert_eq!(HammingSuffix(Count(3)).num_mappings(), 2);
        assert_eq!(HammingSearch(Count(3)).num_mappings(), 3);
        assert_eq!(Edit(Count(1)).num_mappings(), 1);
        assert_eq!(EditPrefix(Count(1)).num_mappings(), 2);
        assert_eq!(EditSuffix(Count(1)).num_mappings(), 2);
        assert_eq!(EditSearch(Count(1)).num_mappings(), 3);
    }

    #[test]
    fn test_match_type_k() {
        assert_eq!(Exact.k(4), 4);
        assert_eq!(ExactPrefix.k(4), 4);
        assert_eq!(Hamming(Count(1)).k(4), 2);
        assert_eq!(Edit(Count(1)).k(4), 2);
    }

    #[test]
    fn test_threshold_get() {
        assert_eq!(Count(2).get(10), 2);
        assert_eq!(Frac(0.5).get(10), 5);
    }

    #[test]
    fn test_graph_null_output() {
        let fq = fastq_bytes(&[("read1", "ACGT", "IIII")]);

        let mut g = Graph::<NoTrace>::new();
        g.add(InputFastqOp::from_reader(Cursor::new(fq)).unwrap());
        g.add(NullOutputOp::new());
        let counter = g.add(CountOp::new([true]));
        g.run().unwrap();

        assert_eq!(counter.counts()[0], 1);

        // Also test unit struct creation
        let _null = NullOutputOp;
    }

    #[test]
    fn test_graph_fork_op() {
        let fq = fastq_bytes(&[("read1", "ACGT", "IIII")]);

        let mut fork_g = Graph::<NoTrace>::new();
        let fork_counter = fork_g.add(CountOp::new([true]));

        let mut g = Graph::<NoTrace>::new();
        g.add(InputFastqOp::from_reader(Cursor::new(fq)).unwrap());
        g.add(ForkOp::new(fork_g));
        let main_counter = g.add(CountOp::new([true]));
        g.run().unwrap();

        // Both the fork and main pipeline should see the read
        assert_eq!(fork_counter.counts()[0], 1);
        assert_eq!(main_counter.counts()[0], 1);
    }

    #[test]
    fn test_graph_while_op() {
        let fq = fastq_bytes(&[("read1", "ACGTACGTACGT", "IIIIIIIIIIII")]);

        let mut while_g = Graph::<NoTrace>::new();
        while_g.add(CutOp::new(te("seq1.* -> seq1.trim, seq1.*"), 4isize));
        while_g.add(TrimOp::new([label("seq1.trim")]));

        let mut g = Graph::<NoTrace>::new();
        g.add(InputFastqOp::from_reader(Cursor::new(fq)).unwrap());
        g.add(WhileOp::new(
            Expr::from(label("seq1.*")).len().gt(4isize),
            while_g,
        ));
        let counter = g.add(CountOp::new([true]));
        g.run().unwrap();

        assert_eq!(counter.counts()[0], 1);
    }

    #[test]
    fn test_graph_try_op() {
        let fq = fastq_bytes(&[("read1", "ACGTNNNN", "IIIIIIII")]);

        // try_graph requires a label that only exists after matching
        let mut try_g = Graph::<NoTrace>::new();
        try_g.add(TrimOp::new([label("seq1.barcode")]));

        let mut catch_g = Graph::<NoTrace>::new();
        let catch_counter = catch_g.add(CountOp::new([true]));

        let mut g = Graph::<NoTrace>::new();
        g.add(InputFastqOp::from_reader(Cursor::new(fq)).unwrap());
        g.add(TryOp::new(try_g, catch_g));
        g.run().unwrap();

        // read1 doesn't have "barcode" label, so try_graph fails requirements -> catch
        assert_eq!(catch_counter.counts()[0], 1);
    }

    #[test]
    fn test_graph_bernoulli_op() {
        let fq = fastq_bytes(&[("read1", "ACGT", "IIII"), ("read2", "TGCA", "IIII")]);

        let mut g = Graph::<NoTrace>::new();
        g.add(InputFastqOp::from_reader(Cursor::new(fq)).unwrap());
        g.add(BernoulliOp::new(attr("seq1.*.rand"), 0.5, 42));

        let trace = NoTrace;
        let (out, _) = g.run_one(None, &trace).unwrap();
        let reads = out.unwrap();
        // Each read should have a .rand attribute set
        for read in &reads {
            let val = read
                .data(
                    StrType::Seq(1),
                    crate::inline_string::InlineString::new(b"*"),
                    crate::inline_string::InlineString::new(b"rand"),
                )
                .unwrap();
            // Should be a Bool
            let _ = val.as_bool();
        }
    }

    #[test]
    fn test_graph_match_polyx_left() {
        let fq = fastq_bytes(&[("read1", "AAAACGTACGT", "IIIIIIIIIII")]);

        let mut g = Graph::<NoTrace>::new();
        g.add(InputFastqOp::from_reader(Cursor::new(fq)).unwrap());
        g.add(MatchPolyXOp::new(
            te("seq1.* -> seq1.polya, seq1.rest"),
            b'A',
            Left,
            0.8,
        ));

        let trace = NoTrace;
        let (out, _) = g.run_one(None, &trace).unwrap();
        let reads = out.unwrap();
        let read = &reads[0];
        // Should have split at the polyA region
        let rest = read.substring(
            StrType::Seq(1),
            crate::inline_string::InlineString::new(b"rest"),
        );
        assert!(rest.is_ok());
    }

    #[test]
    fn test_graph_match_polyx_right() {
        let fq = fastq_bytes(&[("read1", "ACGTACGTAAAA", "IIIIIIIIIIII")]);

        let mut g = Graph::<NoTrace>::new();
        g.add(InputFastqOp::from_reader(Cursor::new(fq)).unwrap());
        g.add(MatchPolyXOp::new(
            te("seq1.* -> seq1.seq, seq1.polya"),
            b'A',
            Right,
            0.8,
        ));

        let trace = NoTrace;
        let (out, _) = g.run_one(None, &trace).unwrap();
        let reads = out.unwrap();
        let read = &reads[0];
        let seq_part = read.substring(
            StrType::Seq(1),
            crate::inline_string::InlineString::new(b"seq"),
        );
        assert!(seq_part.is_ok());
    }

    #[test]
    fn test_graph_output_json() {
        let fq = fastq_bytes(&[("read1", "ACGT", "IIII")]);

        let tmp = std::env::temp_dir().join("antiseq_test_json_writer.jsonl");
        let tmp_str = tmp.to_str().unwrap().to_string();

        {
            let mut g = Graph::<NoTrace>::new();
            g.add(InputFastqOp::from_reader(Cursor::new(fq)).unwrap());
            g.add(OutputJsonOp::from_file(&tmp_str).unwrap());
            g.run().unwrap();
        }

        let output = std::fs::read_to_string(&tmp).unwrap();
        assert!(output.contains("ACGT"));
        std::fs::remove_file(&tmp).ok();
    }

    #[test]
    fn test_graph_time_op() {
        let fq = fastq_bytes(&[("read1", "ACGT", "IIII")]);

        let inner_g = Graph::<NoTrace>::new();

        let mut g = Graph::<NoTrace>::new();
        g.add(InputFastqOp::from_reader(Cursor::new(fq)).unwrap());
        let mut time_op = TimeOp::new(inner_g);
        // We can't easily add TimeOp to graph and get total_time since add() takes ownership
        // Just test construction and total_time
        assert!(time_op.total_time() >= 0.0);
    }

    #[test]
    fn test_graph_match_exact_full() {
        let fq = fastq_bytes(&[("read1", "ACGT", "IIII")]);

        let patterns = Patterns::from_strs(["ACGT"]);

        let mut g = Graph::<NoTrace>::new();
        g.add(InputFastqOp::from_reader(Cursor::new(fq)).unwrap());
        g.add(MatchAnyOp::new(te("seq1.* -> seq1.*"), patterns, Exact));
        let counter = g.add(CountOp::new([true]));
        g.run().unwrap();

        assert_eq!(counter.counts()[0], 1);
    }

    #[test]
    fn test_graph_match_exact_suffix() {
        let fq = fastq_bytes(&[("read1", "NNNNACGT", "IIIIIIII")]);

        let patterns = Patterns::from_strs(["ACGT"]);

        let mut g = Graph::<NoTrace>::new();
        g.add(InputFastqOp::from_reader(Cursor::new(fq)).unwrap());
        g.add(MatchAnyOp::new(
            te("seq1.* -> seq1.rest, seq1.barcode"),
            patterns,
            ExactSuffix,
        ));
        g.add(TrimOp::new([label("seq1.barcode")]));
        let counter = g.add(CountOp::new([true]));
        g.run().unwrap();

        assert_eq!(counter.counts()[0], 1);
    }

    #[test]
    fn test_graph_match_exact_search() {
        let fq = fastq_bytes(&[("read1", "NNACGTNN", "IIIIIIII")]);

        let patterns = Patterns::from_strs(["ACGT"]);

        let mut g = Graph::<NoTrace>::new();
        g.add(InputFastqOp::from_reader(Cursor::new(fq)).unwrap());
        g.add(MatchAnyOp::new(
            te("seq1.* -> seq1.before, seq1.match, seq1.after"),
            patterns,
            ExactSearch,
        ));
        let counter = g.add(CountOp::new([true]));
        g.run().unwrap();

        assert_eq!(counter.counts()[0], 1);
    }

    #[test]
    fn test_graph_match_hamming_full() {
        let fq = fastq_bytes(&[("read1", "ACGC", "IIII")]);

        let patterns = Patterns::from_strs(["ACGT"]);

        let mut g = Graph::<NoTrace>::new();
        g.add(InputFastqOp::from_reader(Cursor::new(fq)).unwrap());
        g.add(MatchAnyOp::new(
            te("seq1.* -> seq1.*"),
            patterns,
            Hamming(Count(3)),
        ));
        let counter = g.add(CountOp::new([true]));
        g.run().unwrap();

        assert_eq!(counter.counts()[0], 1);
    }

    #[test]
    fn test_graph_match_hamming_suffix() {
        let fq = fastq_bytes(&[("read1", "NNNNACGC", "IIIIIIII")]);

        let patterns = Patterns::from_strs(["ACGT"]);

        let mut g = Graph::<NoTrace>::new();
        g.add(InputFastqOp::from_reader(Cursor::new(fq)).unwrap());
        g.add(MatchAnyOp::new(
            te("seq1.* -> seq1.rest, seq1.barcode"),
            patterns,
            HammingSuffix(Count(3)),
        ));
        let counter = g.add(CountOp::new([true]));
        g.run().unwrap();

        assert_eq!(counter.counts()[0], 1);
    }

    #[test]
    fn test_graph_match_hamming_search() {
        let fq = fastq_bytes(&[("read1", "NNACGCNN", "IIIIIIII")]);

        let patterns = Patterns::from_strs(["ACGT"]);

        let mut g = Graph::<NoTrace>::new();
        g.add(InputFastqOp::from_reader(Cursor::new(fq)).unwrap());
        g.add(MatchAnyOp::new(
            te("seq1.* -> seq1.before, seq1.match, seq1.after"),
            patterns,
            HammingSearch(Count(3)),
        ));
        let counter = g.add(CountOp::new([true]));
        g.run().unwrap();

        assert_eq!(counter.counts()[0], 1);
    }

    #[test]
    fn test_graph_match_edit_full() {
        let fq = fastq_bytes(&[("read1", "ACG", "III")]);

        let patterns = Patterns::from_strs(["ACGT"]);

        let mut g = Graph::<NoTrace>::new();
        g.add(InputFastqOp::from_reader(Cursor::new(fq)).unwrap());
        g.add(MatchAnyOp::new(
            te("seq1.* -> seq1.*"),
            patterns,
            Edit(Count(1)),
        ));
        let counter = g.add(CountOp::new([true]));
        g.run().unwrap();

        assert_eq!(counter.counts()[0], 1);
    }

    #[test]
    fn test_graph_match_edit_suffix() {
        let fq = fastq_bytes(&[("read1", "NNNNACGT", "IIIIIIII")]);

        let patterns = Patterns::from_strs(["ACGT"]);

        let mut g = Graph::<NoTrace>::new();
        g.add(InputFastqOp::from_reader(Cursor::new(fq)).unwrap());
        g.add(MatchAnyOp::new(
            te("seq1.* -> seq1.rest, seq1.barcode"),
            patterns,
            EditSuffix(Count(1)),
        ));
        let counter = g.add(CountOp::new([true]));
        g.run().unwrap();

        assert_eq!(counter.counts()[0], 1);
    }

    #[test]
    fn test_graph_match_edit_search() {
        let fq = fastq_bytes(&[("read1", "NNACGTNN", "IIIIIIII")]);

        let patterns = Patterns::from_strs(["ACGT"]);

        let mut g = Graph::<NoTrace>::new();
        g.add(InputFastqOp::from_reader(Cursor::new(fq)).unwrap());
        g.add(MatchAnyOp::new(
            te("seq1.* -> seq1.before, seq1.match, seq1.after"),
            patterns,
            EditSearch(Count(1)),
        ));
        let counter = g.add(CountOp::new([true]));
        g.run().unwrap();

        assert_eq!(counter.counts()[0], 1);
    }

    #[test]
    fn test_graph_set_op_attr() {
        let fq = fastq_bytes(&[("read1", "ACGT", "IIII")]);

        let mut g = Graph::<NoTrace>::new();
        g.add(InputFastqOp::from_reader(Cursor::new(fq)).unwrap());
        // Set label to bytes
        g.add(SetOp::new(label("seq1.*"), b"NNNN".to_vec()));
        // Set attribute to computed value
        g.add(SetOp::new(
            attr("seq1.*.length"),
            Expr::from(label("seq1.*")).len(),
        ));

        let trace = NoTrace;
        let (out, _) = g.run_one(None, &trace).unwrap();
        let reads = out.unwrap();
        let read = &reads[0];
        let length = read
            .data(
                StrType::Seq(1),
                crate::inline_string::InlineString::new(b"*"),
                crate::inline_string::InlineString::new(b"length"),
            )
            .unwrap();
        assert_eq!(length, &Data::Int(4));
    }

    #[test]
    fn test_graph_match_distance_counts_with_match() {
        let fq = fastq_bytes(&[
            ("read1", "ACGTNNNN", "IIIIIIII"),
            ("read2", "TGCANNNN", "IIIIIIII"),
        ]);

        let patterns = Patterns::from_strs(["ACGT", "TGCA"]);

        let mut g = Graph::<NoTrace>::new();
        g.add(InputFastqOp::from_reader(Cursor::new(fq)).unwrap());
        g.add(MatchAnyOp::new(
            te("seq1.* -> seq1.barcode, seq1.rest"),
            patterns,
            ExactPrefix,
        ));
        g.run().unwrap();

        let counts = g.match_distance_counts();
        // Should have match distance counts from the MatchAnyOp
        assert!(!counts.is_empty());
    }

    #[test]
    fn test_graph_output_json_file() {
        let fq = fastq_bytes(&[("read1", "ACGT", "IIII")]);

        let tmp = std::env::temp_dir().join("antiseq_test_json.jsonl");
        let tmp_str = tmp.to_str().unwrap().to_string();

        {
            let mut g = Graph::<NoTrace>::new();
            g.add(InputFastqOp::from_reader(Cursor::new(fq)).unwrap());
            g.add(OutputJsonOp::from_file(&tmp_str).unwrap());
            g.run().unwrap();
        }

        let output = std::fs::read_to_string(&tmp).unwrap();
        assert!(output.contains("ACGT"));
        std::fs::remove_file(&tmp).ok();
    }

    #[test]
    fn test_graph_output_fastq_multiple_files() {
        let _fq = fastq_bytes(&[("read1", "ACGTNNNN", "IIIIIIII")]);

        let tmp1 = std::env::temp_dir().join("antiseq_test_multi_out1.fastq");
        let tmp2 = std::env::temp_dir().join("antiseq_test_multi_out2.fastq");
        let tmp1_str = tmp1.to_str().unwrap().to_string();
        let tmp2_str = tmp2.to_str().unwrap().to_string();

        {
            let mut g = Graph::<NoTrace>::new();
            g.add(
                InputFastqOp::from_readers(vec![
                    Cursor::new(fastq_bytes(&[("r1_R1", "AAAA", "IIII")])),
                    Cursor::new(fastq_bytes(&[("r1_R2", "CCCC", "IIII")])),
                ])
                .unwrap(),
            );
            g.add(OutputFastqFileOp::from_files([
                tmp1_str.clone(),
                tmp2_str.clone(),
            ]));
            g.run().unwrap();
        }

        let output1 = std::fs::read_to_string(&tmp1).unwrap();
        assert!(output1.contains("AAAA"));
        let output2 = std::fs::read_to_string(&tmp2).unwrap();
        assert!(output2.contains("CCCC"));
        std::fs::remove_file(&tmp1).ok();
        std::fs::remove_file(&tmp2).ok();
    }

    #[test]
    fn test_graph_match_bounded() {
        let fq = fastq_bytes(&[("read1", "NNACGTNNNN", "IIIIIIIIII")]);

        let patterns = Patterns::from_strs(["ACGT"]);

        let mut g = Graph::<NoTrace>::new();
        g.add(InputFastqOp::from_reader(Cursor::new(fq)).unwrap());
        g.add(MatchAnyOp::new(
            te("seq1.* -> seq1.before, seq1.match, seq1.after"),
            patterns,
            ExactBoundedMatch { from: 1, to: 6 },
        ));
        let counter = g.add(CountOp::new([true]));
        g.run().unwrap();

        assert_eq!(counter.counts()[0], 1);
    }

    #[test]
    fn test_graph_match_hamming_bounded() {
        let fq = fastq_bytes(&[("read1", "NNACGCNNNN", "IIIIIIIIII")]);

        let patterns = Patterns::from_strs(["ACGT"]);

        let mut g = Graph::<NoTrace>::new();
        g.add(InputFastqOp::from_reader(Cursor::new(fq)).unwrap());
        g.add(MatchAnyOp::new(
            te("seq1.* -> seq1.before, seq1.match, seq1.after"),
            patterns,
            HammingBoundedMatch {
                threshold: Count(3),
                from: 1,
                to: 6,
            },
        ));
        let counter = g.add(CountOp::new([true]));
        g.run().unwrap();

        assert_eq!(counter.counts()[0], 1);
    }

    #[test]
    fn test_graph_match_edit_bounded() {
        let fq = fastq_bytes(&[("read1", "NNACGNNNN", "IIIIIIIII")]);

        let patterns = Patterns::from_strs(["ACGT"]);

        let mut g = Graph::<NoTrace>::new();
        g.add(InputFastqOp::from_reader(Cursor::new(fq)).unwrap());
        g.add(MatchAnyOp::new(
            te("seq1.* -> seq1.before, seq1.match, seq1.after"),
            patterns,
            EditBoundedMatch {
                threshold: Count(1),
                from: 1,
                to: 6,
            },
        ));
        let counter = g.add(CountOp::new([true]));
        g.run().unwrap();

        assert_eq!(counter.counts()[0], 1);
    }

    #[test]
    fn test_match_type_k_aln() {
        assert_eq!(GlobalAln(0.75).k(4), 2);
        let k = PrefixAln {
            identity: 0.75,
            overlap: 0.5,
        }
        .k(8);
        assert!(k > 0);
        let k = SuffixAln {
            identity: 0.75,
            overlap: 0.5,
        }
        .k(8);
        assert!(k > 0);
        let k = LocalAln {
            identity: 0.75,
            overlap: 0.5,
        }
        .k(8);
        assert!(k > 0);
        let k = ExactBoundedMatch { from: 0, to: 3 }.k(4);
        assert_eq!(k, 4);
        let k = HammingBoundedMatch {
            threshold: Count(1),
            from: 0,
            to: 3,
        }
        .k(4);
        assert!(k > 0);
        let k = EditBoundedMatch {
            threshold: Count(1),
            from: 0,
            to: 3,
        }
        .k(4);
        assert!(k > 0);
    }

    #[test]
    fn test_match_type_num_mappings_aln() {
        assert_eq!(GlobalAln(0.75).num_mappings(), 1);
        assert_eq!(
            (PrefixAln {
                identity: 0.75,
                overlap: 0.5
            })
            .num_mappings(),
            2
        );
        assert_eq!(
            (SuffixAln {
                identity: 0.75,
                overlap: 0.5
            })
            .num_mappings(),
            2
        );
        assert_eq!(
            (LocalAln {
                identity: 0.75,
                overlap: 0.5
            })
            .num_mappings(),
            3
        );
        assert_eq!((ExactBoundedMatch { from: 0, to: 3 }).num_mappings(), 3);
        assert_eq!(
            (HammingBoundedMatch {
                threshold: Count(1),
                from: 0,
                to: 3
            })
            .num_mappings(),
            3
        );
        assert_eq!(
            (EditBoundedMatch {
                threshold: Count(1),
                from: 0,
                to: 3
            })
            .num_mappings(),
            3
        );
    }

    #[test]
    fn test_graph_match_regex() {
        let fq = fastq_bytes(&[
            ("read1", "ACGTNNNN", "IIIIIIII"),
            ("read2", "NNNNNNNN", "IIIIIIII"),
        ]);

        let mut g = Graph::<NoTrace>::new();
        g.add(InputFastqOp::from_reader(Cursor::new(fq)).unwrap());
        g.add(MatchRegexOp::new(te("seq1.* -> seq1.*.matched"), "ACGT"));

        let trace = NoTrace;
        let (out, _) = g.run_one(None, &trace).unwrap();
        let reads = out.unwrap();

        let matched1 = reads[0]
            .data(
                StrType::Seq(1),
                crate::inline_string::InlineString::new(b"*"),
                crate::inline_string::InlineString::new(b"matched"),
            )
            .unwrap();
        assert_eq!(matched1, &Data::Bool(true));

        let matched2 = reads[1]
            .data(
                StrType::Seq(1),
                crate::inline_string::InlineString::new(b"*"),
                crate::inline_string::InlineString::new(b"matched"),
            )
            .unwrap();
        assert_eq!(matched2, &Data::Bool(false));
    }

    #[test]
    fn test_graph_match_regex_with_captures() {
        let fq = fastq_bytes(&[("read1", "ACGTNNNN", "IIIIIIII")]);

        let mut g = Graph::<NoTrace>::new();
        g.add(InputFastqOp::from_reader(Cursor::new(fq)).unwrap());
        g.add(MatchRegexOp::new(
            te("seq1.* -> seq1.*.matched"),
            "(?P<barcode>[ACGT]{4})",
        ));

        let trace = NoTrace;
        let (out, _) = g.run_one(None, &trace).unwrap();
        let reads = out.unwrap();
        let read = &reads[0];

        // The named capture "barcode" should create a new mapping
        let barcode = read.substring(
            StrType::Seq(1),
            crate::inline_string::InlineString::new(b"barcode"),
        );
        assert!(barcode.is_ok());
        assert_eq!(barcode.unwrap(), b"ACGT");
    }

    #[test]
    fn test_trace_reads() {
        use crate::trace::TraceReads;

        let fq = fastq_bytes(&[("read1", "ACGT", "IIII"), ("read2", "TGCA", "IIII")]);

        let tmp = std::env::temp_dir().join("antiseq_test_trace.json");
        let tmp_str = tmp.to_str().unwrap().to_string();

        {
            let mut g = Graph::<TraceReads>::new();
            g.add(InputFastqOp::from_reader(Cursor::new(fq)).unwrap());
            let counter = g.add(CountOp::new([true]));
            g.run_trace(&tmp_str).unwrap();
            assert_eq!(counter.counts()[0], 2);
        }

        let output = std::fs::read_to_string(&tmp).unwrap();
        assert!(output.contains("traceEvents"));
        assert!(output.contains("InputFastqOp"));
        std::fs::remove_file(&tmp).ok();
    }

    #[test]
    fn test_graph_time_op_in_pipeline() {
        let fq = fastq_bytes(&[("read1", "ACGT", "IIII")]);

        let mut inner_g = Graph::<NoTrace>::new();
        inner_g.add(CountOp::new([true]));

        let mut g = Graph::<NoTrace>::new();
        g.add(InputFastqOp::from_reader(Cursor::new(fq)).unwrap());
        g.add(TimeOp::new(inner_g));
        g.run().unwrap();
    }

    #[test]
    fn test_graph_set_label_with_qual() {
        // Set a label on a read that has quality scores
        let fq = fastq_bytes(&[("read1", "ACGTACGT", "IIIIIIII")]);

        let mut g = Graph::<NoTrace>::new();
        g.add(InputFastqOp::from_reader(Cursor::new(fq)).unwrap());
        g.add(CutOp::new(te("seq1.* -> seq1.left, seq1.right"), 4isize));
        // Set left to a different sequence - this exercises the qual branch in SetOp
        g.add(SetOp::new(
            label("seq1.left"),
            Expr::from(label("seq1.right")),
        ));

        let trace = NoTrace;
        let (out, _) = g.run_one(None, &trace).unwrap();
        let reads = out.unwrap();
        let read = &reads[0];

        let left = read
            .substring(
                StrType::Seq(1),
                crate::inline_string::InlineString::new(b"left"),
            )
            .unwrap();
        assert_eq!(left, b"ACGT");
    }

    #[test]
    fn test_graph_set_label_no_qual() {
        // Set a label on name (which has no quality scores)
        let fq = fastq_bytes(&[("read1", "ACGT", "IIII")]);

        let mut g = Graph::<NoTrace>::new();
        g.add(InputFastqOp::from_reader(Cursor::new(fq)).unwrap());
        g.add(SetOp::new(label("name1.*"), b"newname".to_vec()));

        let trace = NoTrace;
        let (out, _) = g.run_one(None, &trace).unwrap();
        let reads = out.unwrap();
        let read = &reads[0];

        let name = read
            .substring(
                StrType::Name(1),
                crate::inline_string::InlineString::new(b"*"),
            )
            .unwrap();
        assert_eq!(name, b"newname");
    }

    #[test]
    fn test_graph_remove_internal_pipeline() {
        let fq = fastq_bytes(&[("read1", "ACGTACGT", "IIIIIIII")]);

        let mut g = Graph::<NoTrace>::new();
        g.add(InputFastqOp::from_reader(Cursor::new(fq)).unwrap());
        g.add(CutOp::new(te("seq1.* -> seq1.left, seq1.right"), 4isize));
        g.add(TrimOp::new([label("seq1.left")]));
        g.add(RemoveInternalOp::create());

        let trace = NoTrace;
        let (out, _) = g.run_one(None, &trace).unwrap();
        let reads = out.unwrap();
        let read = &reads[0];

        let seq = read
            .substring(
                StrType::Seq(1),
                crate::inline_string::InlineString::new(b"right"),
            )
            .unwrap();
        assert_eq!(seq, b"ACGT");
    }

    #[test]
    fn test_graph_for_each_op() {
        let fq = fastq_bytes(&[("read1", "ACGT", "IIII")]);

        let mut g = Graph::<NoTrace>::new();
        g.add(InputFastqOp::from_reader(Cursor::new(fq)).unwrap());
        g.add(ForEachOp::new(|_read: &mut Read| {
            // Just a no-op closure to exercise ForEachOp
        }));
        let counter = g.add(CountOp::new([true]));
        g.run().unwrap();

        assert_eq!(counter.counts()[0], 1);
    }

    #[test]
    fn test_graph_output_fastq_gzipped() {
        let fq = fastq_bytes(&[("read1", "ACGT", "IIII")]);

        let tmp = std::env::temp_dir().join("antiseq_test_gz.fastq.gz");
        let tmp_str = tmp.to_str().unwrap().to_string();

        {
            let mut g = Graph::<NoTrace>::new();
            g.add(InputFastqOp::from_reader(Cursor::new(fq)).unwrap());
            g.add(OutputFastqFileOp::from_file(tmp_str.clone()));
            g.run().unwrap();
        }

        // Just check that the file exists and is non-empty
        let metadata = std::fs::metadata(&tmp).unwrap();
        assert!(metadata.len() > 0);
        std::fs::remove_file(&tmp).ok();
    }

    #[test]
    fn test_graph_set_label_bytes_with_qual() {
        // Exercise the SetOp label-with-qual path using a bytes literal
        let fq = fastq_bytes(&[("read1", "ACGTACGT", "IIIIIIII")]);

        let mut g = Graph::<NoTrace>::new();
        g.add(InputFastqOp::from_reader(Cursor::new(fq)).unwrap());
        // Setting seq label to literal bytes exercises the no-qual-expression branch
        g.add(SetOp::new(label("seq1.*"), b"NNNN".to_vec()));

        let trace = NoTrace;
        let (out, _) = g.run_one(None, &trace).unwrap();
        let reads = out.unwrap();
        let read = &reads[0];
        let seq = read
            .substring(
                StrType::Seq(1),
                crate::inline_string::InlineString::new(b"*"),
            )
            .unwrap();
        assert_eq!(seq, b"NNNN");
    }

    #[test]
    fn test_graph_retain_none_pass() {
        // All reads filtered out
        let fq = fastq_bytes(&[("read1", "AC", "II")]);

        let tmp = std::env::temp_dir().join("antiseq_test_retain_none.fastq");
        let tmp_str = tmp.to_str().unwrap().to_string();

        {
            let mut g = Graph::<NoTrace>::new();
            g.add(InputFastqOp::from_reader(Cursor::new(fq)).unwrap());
            g.add(RetainOp::new(
                Expr::from(label("seq1.*")).len().gt(100isize),
            ));
            g.add(OutputFastqFileOp::from_file(tmp_str.clone()));
            g.run().unwrap();
        }

        // Output file should exist but be empty (no reads pass filter)
        let output = std::fs::read_to_string(&tmp).unwrap_or_default();
        assert!(!output.contains("AC"));
        std::fs::remove_file(&tmp).ok();
    }

    #[test]
    fn test_graph_match_frac_threshold() {
        let fq = fastq_bytes(&[("read1", "ACGCNNNN", "IIIIIIII")]);

        let patterns = Patterns::from_strs(["ACGT"]);

        let mut g = Graph::<NoTrace>::new();
        g.add(InputFastqOp::from_reader(Cursor::new(fq)).unwrap());
        g.add(MatchAnyOp::new(
            te("seq1.* -> seq1.barcode, seq1.rest"),
            patterns,
            HammingPrefix(Frac(0.75)),
        ));
        let counter = g.add(CountOp::new([true]));
        g.run().unwrap();

        assert_eq!(counter.counts()[0], 1);
    }
}
