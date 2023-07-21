use antisequence::*;

fn main() {
    iter_fastq1("example_data/filter_l.fastq", 256)
        .unwrap_or_else(|e| panic!("{e}"))
        .dbg(sel!())
        .filter(
            sel!(),
            tr!(seq1.* -> seq1.*._f),
            String::from("example_data/filter.txt"),
            2,
        )
        .dbg(sel!())
        .run()
        .unwrap_or_else(|e| panic!("{e}"));
}
// GAGGTATAGTCTTG
// GAGGTATGGTCTTG
// GAGGTAT
//  AGGTATA
//   GGTATAG
