use antisequence::expr::*;
use antisequence::node::*;
use antisequence::*;

fn main() {
    // Demo single-cell sequencing protocol:
    // R1: bc[9-11] CAGAGC umi[8] bc[10]
    // R2: insert adapter

    let mut g = Graph::new();
    g.add(
        InputFastq1Node::new_interleaved("./example_data/single_cell.fastq")
            .unwrap_or_else(|e| panic!("{e}")),
    );
    g.add(DbgNode::new());
    g.add(MatchAnyNode::new(
        tr!(seq2.* -> _, seq2.adapter),
        Patterns::from_exprs([Expr::from("ATATATATAT"), Expr::from("CGCGCGCGCG")]),
        SuffixAln {
            identity: 0.7,
            overlap: 0.4,
        },
    ));
    g.add(DbgNode::new());
    g.add(TrimNode::new([label("seq2.adapter")]));
    g.add(MatchAnyNode::new(
        tr!(seq1.* -> seq1.bc1, _, seq1.after_anchor),
        Patterns::from_exprs([Expr::from("CAGAGC")]),
        HammingSearch(Frac(0.8)),
    ));
    g.add(CutNode::new(
        tr!(seq1.after_anchor -> seq1.umi, seq1.after_umi),
        LeftEnd(8),
    ));
    g.add(CutNode::new(
        tr!(seq1.after_umi -> seq1.bc2, _),
        LeftEnd(10),
    ));
    g.add(RetainNode::new(
        Expr::from(label_exists("seq1.bc1")).and(Expr::from(label_exists("seq1.bc2"))),
    ));
    g.add(RetainNode::new(
        Expr::from(label("seq1.brc1"))
            .in_bounds(Expr::from(9 as usize)..=Expr::from(11 as usize))
            .and(
                Expr::from(label("seq1.brc2"))
                    .len()
                    .eq(Expr::from(10 as usize)),
            ),
    ));
    g.add(SetNode::new(
        label("name1.*"),
        concat_all([
            Expr::from(label("name1.*")),
            Expr::from("_"),
            Expr::from(label("seq1.umi")),
            Expr::from("_"),
            Expr::from(label("seq1.bc1")),
            Expr::from(label("seq1.bc2")),
        ]),
    ));
    g.add(SetNode::new(label("seq1.*"), Expr::from(label("seq2.*"))));
    g.add(DbgNode::new());
    g.add(OutputFastqNode::new1("example_output/single_cell.fastq"));
    g.run().unwrap_or_else(|e| panic!("{e}"));
}
