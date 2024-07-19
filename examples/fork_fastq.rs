use antisequence::expr::*;
use antisequence::node::*;
use antisequence::*;

fn main() {
    let mut right = Graph::new();
    right.add(InputFastq1Node::new("example_data/simple.fastq").unwrap_or_else(|e| panic!("{e}")));
    right.add(CutNode::new(tr!(seq1.* -> seq1.a, seq1.b), LeftEnd(3)));
    right.add(TrimNode::new([label("seq1.a")]));
    right.add(DbgNode::new());

    let mut left = Graph::new();
    left.add(SetNode::new(
        label("name1.*"),
        concat_all([Expr::from(label("name1.*")), Expr::from(label("seq1.a"))]),
    ));
    left.add(DbgNode::new());
    right.add(ForkNode::new(left));

    right.run().unwrap_or_else(|e| panic!("{e}"));
}
