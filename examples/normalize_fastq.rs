use antisequence::expr::*;
use antisequence::node::*;
use antisequence::*;

fn main() {
    let fastq = b"@read1
TTTTTTGAC
+
012345678
@read1a
TTTTTTGAC
+
012345678
@read2
AAAAAAAGAC
+
0123456789
@read3
CCCCCCCCGAC
+
01234567890";

    let mut g = Graph::new();
    g.add(InputFastq1Node::from_bytes(fastq).unwrap_or_else(|e| panic!("{e}")));
    g.add(MatchAnyNode::new(
        tr!(seq1.* -> seq1.left, seq1.patt),
        Patterns::from_strs(["GAC"]),
        ExactSuffix,
    ));
    g.add(SetNode::new(
        label("seq1.left"),
        Expr::from(label("seq1.left")).normalize(6..=8),
    ));
    g.add(DbgNode::new());
    g.add(OutputFastqNode::new1("example_output/simple.fastq"));
    g.run().unwrap_or_else(|e| panic!("{e}"));
}
