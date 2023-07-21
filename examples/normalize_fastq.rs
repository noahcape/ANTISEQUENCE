use antisequence::*;

fn main() {
    iter_fastq1("example_data/normalize.fastq", 256)
        .unwrap_or_else(|e| panic!("{e}"))
        .match_one(sel!(), tr!(seq1.* -> seq1.l, seq1.anchor, seq1.rest), "CAGAGC", HammingSearch(Frac(1.0)))
        .dbg(sel!(seq1.anchor))
        .norm(sel!(seq1.l), label!(seq1.l), 6..=11)
        .dbg(sel!(seq1.l))
        .run()
        .unwrap_or_else(|e| panic!("{e}"));
}
