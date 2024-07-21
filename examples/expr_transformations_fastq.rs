use antisequence::expr::*;
use antisequence::node::*;
use antisequence::*;

fn main() {
    let fastq = b"@read1/1
AAAAACCATTTTT
+
0123456789012
@read1/2
AAAAACCATTTTT
+
0123456789012";

    let mut g = Graph::new();
    g.add(InputFastq1Node::from_bytes(fastq).unwrap_or_else(|e| panic!("{e}")));

    g.add(CutNode::new(tr!(seq1.* -> seq1.a, seq1.b), LeftEnd(5)));
    g.add(CutNode::new(tr!(seq1.b -> seq1.mid, seq1.bb), LeftEnd(3)));
    g.add(SetNode::new(
        label("seq1.mid"),
        Expr::from(label("seq1.mid")).rev(),
    ));
    g.add(SetNode::new(
        label("seq1.mid"),
        Expr::from(label("seq1.mid")).slice(Expr::from(1)..),
    ));
    g.add(SetNode::new(
        label("seq1.mid"),
        Expr::from(label("seq1.mid")).revcomp(),
    ));
    g.add(DbgNode::new());
    g.add(OutputFastqNode::new1("example_output/simple.fastq"));
    g.run().unwrap_or_else(|e| panic!("{e}"));
}
