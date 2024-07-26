use antisequence::expr::*;
use antisequence::graph::*;
use antisequence::*;

fn main() {
    let mut g = Graph::new();
    g.add(InputFastqOp::from_file("example_data/simple.fastq").unwrap_or_else(|e| panic!("{e}")));
    g.add(CutOp::new(tr!(seq1.* -> seq1.a, seq1.b), LeftEnd(3)));
    g.add(CutOp::new(tr!(seq1.b -> _, seq1.b), RightEnd(4)));
    g.add(DbgOp::new());
    g.add(SetOp::new(label("name1.*"), fmt_expr("{name1.*}_{seq1.a}")));
    g.add(TrimOp::new([label("seq1.a")]));
    g.add(DbgOp::new());
    g.add(OutputFastqFileOp::from_file("example_output/simple.fastq"));
    g.add(OutputJsonOp::from_file("example_output/simple.json").unwrap_or_else(|e| panic!("{e}")));
    g.run().unwrap_or_else(|e| panic!("{e}"));
}
