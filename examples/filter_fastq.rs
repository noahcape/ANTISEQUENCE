use antisequence::*;

fn main() {
    iter_fastq1("example_data/filter_l.fastq", 256)
        .unwrap_or_else(|e| panic!("{e}"))
        .match_one(sel!(), tr!(seq1.* -> _, seq1.anchor, seq1.rest), "CAGAGC", HammingSearch(Frac(0.8)))
        .cut(sel!(seq1.anchor), tr!(seq1.rest -> _, seq1.rest_r), LeftEnd(8))
        .cut(sel!(seq1.rest_r), tr!(seq1.rest_r -> seq1.brc, _), LeftEnd(10))
        .filter(
            sel!(seq1.brc),
            tr!(seq1.brc -> seq1.brc._f),
            String::from("example_data/filter.txt"),
            1,
        )
        .dbg(sel!(seq1.brc))
        .run()
        .unwrap_or_else(|e| panic!("{e}"));
}
