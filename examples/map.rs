use antisequence::*;

fn main() {
    iter_fastq1("example_data/map.fastq", 256)
        .unwrap_or_else(|e| panic!("{e}"))
        .cut(sel!(), tr!(seq1.* -> seq1.brc, seq1.rest), LeftEnd(8))
        .map(
            sel!(),
            tr!(seq1.brc -> seq1.brc.not_mapped),
            String::from("example_data/bc_map.txt"),
            1,
        )
        .dbg(sel!())
        .collect_fastq1(sel!(), "example_output/map.fastq")
        .run()
        .unwrap_or_else(|e| panic!("{e}"));
}
