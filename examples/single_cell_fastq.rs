use antisequence::expr::*;
use antisequence::graph::*;
use antisequence::*;

fn main() {
    // Demo single-cell sequencing protocol:
    // R1: bc[9-11] CAGAGC umi[8] bc[10]
    // R2: insert adapter

    let fastq = b"@read1/1
AAAAAAAAAACAGAGCTTTTTTTTCCCCCCCCCC
+
0123456789012345678901234567890123
@read1/2
AAAATTTTCCCCGGGGAAAACGCGACG
+
012345678901234567890123456
@read2/1
AAAAAAAAAAAAAACAGAGCTTTTTTTTCCCCCCCCCC
+
01234567890123456789012345678901234567
@read2/2
AAAATTTTCCCCGGGGATATAT
+
0123456789012345678901";

    let adapters = ["ATATATATAT", "CGCGCGCGCG"];

    let mut g = Graph::new();
    g.add(
        InputFastqNode::from_interleaved_reader(fastq.as_slice(), 2)
            .unwrap_or_else(|e| panic!("{e}")),
    );

    // trim adapter
    g.add(MatchAnyNode::new(
        tr!(seq2.* -> _, seq2.adapter),
        Patterns::from_strs(adapters),
        SuffixAln {
            identity: 0.7,
            overlap: 0.4,
        },
    ));
    g.add(DbgNode::new());
    g.add(TrimNode::new([label("seq2.adapter")]));

    // match anchor
    g.add(MatchAnyNode::new(
        tr!(seq1.* -> seq1.bc1, _, seq1._after_anchor),
        Patterns::from_strs(["CAGAGC"]),
        HammingSearch(Frac(0.8)),
    ));

    // split the UMI from the rest of the sequence
    g.add(CutNode::new(
        tr!(seq1._after_anchor -> seq1.umi, seq1._after_umi),
        LeftEnd(8),
    ));

    // clip the length of the second barcode
    g.add(CutNode::new(
        tr!(seq1._after_umi -> seq1.bc2, _),
        LeftEnd(10),
    ));
    g.add(DbgNode::new());

    // filter out invalid reads
    g.add(RetainNode::new(
        label_exists("seq1.bc1")
            .and(label_exists("seq1.bc2"))
            .and(Expr::from(label("seq1.bc1")).len().in_bounds(9..=11))
            .and(Expr::from(label("seq1.bc2")).len().eq(10)),
    ));

    // move the UMI and barcodes to the read name
    g.add(SetNode::new(
        label("name1.*"),
        fmt_expr("{name1.*}_{seq1.umi}_{seq1.bc1}{seq1.bc2}"),
    ));
    g.add(SetNode::new(label("seq1.*"), label("seq2.*")));

    g.add(OutputFastqFileNode::from_file(
        "example_output/single_cell.fastq",
    ));
    g.run().unwrap_or_else(|e| panic!("{e}"));
}
