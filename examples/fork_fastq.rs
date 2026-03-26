use antisequence::expr::*;
use antisequence::graph::*;
use antisequence::*;

fn main() {
    let mut g = <Graph>::new();
    g.add(InputFastqOp::from_file("example_data/simple.fastq").unwrap_or_else(|e| panic!("{e}")));
    g.add(CutOp::new(tr!(seq1.* -> seq1.a, seq1.b), 3));

    let mut fork = <Graph>::new();
    fork.add(SetOp::new(label("name1.*"), fmt_expr("{name1.*}_{seq1.a}")));
    fork.add(DbgOp::create());
    g.add(ForkOp::new(fork));

    g.add(TrimOp::new([label("seq1.a")]));
    g.add(DbgOp::create());
    g.run().unwrap_or_else(|e| panic!("{e}"));
}
