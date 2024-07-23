use antisequence::expr::*;
use antisequence::graph::*;
use antisequence::*;

fn main() {
    let mut g = Graph::new();
    g.add(InputFastq1Node::new("example_data/simple.fastq").unwrap_or_else(|e| panic!("{e}")));
    g.add(CutNode::new(tr!(seq1.* -> seq1.a, seq1.b), LeftEnd(3)));
    g.add(CutNode::new(tr!(seq1.b -> _, seq1.b), RightEnd(4)));
    g.add(DbgNode::new());
    g.add(SetNode::new(
        label("name1.*"),
        fmt_expr("{name1.*}_{seq1.a}"),
    ));
    g.add(TrimNode::new([label("seq1.a")]));
    g.add(DbgNode::new());
    g.add(OutputFastqNode::new1("example_output/simple.fastq"));
    g.add(ForEachNode::new(|read| println!("{}", read.to_json())));
    g.run().unwrap_or_else(|e| panic!("{e}"));
}
