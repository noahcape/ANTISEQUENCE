use antisequence::*;

fn main() {
    iter_fastq1("./example_data/match_one.fastq", 256)
        .unwrap_or_else(|e| panic!("{e}"))
        .match_one(
            sel!(),
            tr!(seq1.* -> seq1.barcode_2, seq1.anchor, seq1._r),
            "CAGAGC",
            HammingSearch(Frac(0.83)),
        )
        .dbg(sel!())
        .collect_fastq1(sel!(), "./example_output/match_one.fastq")
        .run()
        .unwrap_or_else(|e| panic!("{e}"));
}