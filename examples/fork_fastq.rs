use antisequence::expr::*;
use antisequence::graph::*;
use antisequence::*;

fn main() {
    let mut g = Graph::new();
    g.add(InputFastq1Node::new("example_data/simple.fastq").unwrap_or_else(|e| panic!("{e}")));
    g.add(CutNode::new(tr!(seq1.* -> seq1.a, seq1.b), LeftEnd(3)));

    let mut fork = Graph::new();
    fork.add(SetNode::new(
        label("name1.*"),
        fmt_expr("{name1.*}_{seq1.a}"),
    ));
    fork.add(DbgNode::new());
    g.add(ForkNode::new(fork));

    g.add(TrimNode::new([label("seq1.a")]));
    g.add(DbgNode::new());
    g.run().unwrap_or_else(|e| panic!("{e}"));
}
