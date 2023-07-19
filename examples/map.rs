use antisequence::*;

fn main() {
    iter_fastq1("example_data/map.fastq", 256)
        .unwrap_or_else(|e| panic!("{e}"))
        .dbg(sel!())
        .map(
            sel!(),
            tr!(seq1.* -> seq1.*.mapped),
            String::from("example_data/bc_map.txt"),
            0,
        )
        .pad(sel!(!seq1.*.mapped), [label!(seq1.*)], RightEnd(10), b'0')
        .dbg(sel!())
        .collect_fastq1(sel!(), "example_output/map.fastq")
        .run()
        .unwrap_or_else(|e| panic!("{e}"));
}
