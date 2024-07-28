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
//! ```
//! @read6
//! AATTCCGGAATTCCCAAAAG
//! +
//! 01234567890123456789
//! ```
//! The first, second, and fourth lines are the name, sequence, and quality scores, respectively.
//!
//! ANTISEQUENCE stores that record as an internal [`Read`] data structure:
//! ```
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

pub mod errors;
pub mod expr;
pub mod graph;
mod patterns;
mod read;

mod inline_string;
mod parse_utils;

// commonly used functions and types

pub use crate::patterns::*;
pub use crate::read::*;
