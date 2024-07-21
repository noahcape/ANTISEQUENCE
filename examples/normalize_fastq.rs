use antisequence::expr::*;
use antisequence::node::*;
use antisequence::*;

fn main() {
    let fastq = b"@read1
TTTTTT
+
012345
@read1a
TTTTTT
+
012345
@read2
AAAAAAA
+
0123456
@read3
CCCCCCCC
+
01234567";

    let mut g = Graph::new();
    g.add(InputFastq1Node::from_bytes(fastq).unwrap_or_else(|e| panic!("{e}")));
    g.add(SetNode::new(
        label("seq1.*"),
        Expr::from(label("seq1.*")).normalize(6..=8),
    ));
    g.add(DbgNode::new());
    g.add(OutputFastqNode::new1("example_output/simple.fastq"));
    g.run().unwrap_or_else(|e| panic!("{e}"));
}
