use antisequence::expr::*;
use antisequence::graph::*;
use antisequence::*;

fn main() {
    let mut g = Graph::new();
    g.add(InputFastqNode::from_file("example_data/match.fastq").unwrap_or_else(|e| panic!("{e}")));
    let patterns = Patterns::from_strs(["AAAA", "TTTT"]);
    g.add(MatchAnyNode::new(
        tr!(seq1.* -> seq1.template, seq1.adapter),
        patterns,
        SuffixAln {
            identity: 0.75,
            overlap: 0.5,
        },
    ));
    g.add(DbgNode::new());
    g.add(TrimNode::new([label("seq1.adapter")]));
    g.add(DbgNode::new());
    g.add(OutputFastqFileNode::from_file("example_output/match.fastq"));
    g.run().unwrap_or_else(|e| panic!("{e}"));
}
