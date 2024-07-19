use antisequence::*;

fn main() {
    let file_1 = "./example_data/SRR6750041_1.rand1000.fastq";
    let file_2 = "./example_data/SRR6750041_2.rand1000.fastq";
    let map_file = "./example_data/split_seq_barcodes.txt";

    iter_fastq2(file_1, file_2, 256)
        .unwrap_or_else(|e| panic!("{e}"))
        .cut(sel!(seq2.*), tr!(seq2.* -> seq2._l, seq2.brc1), RightEnd(8))
        .map(
            sel!(seq2.brc1),
            tr!(seq2.brc1 -> seq2.brc1.not_mapped),
            map_file,
            1,
        )
        .cut(
            sel!(seq2._l),
            tr!(seq2._l -> seq2._l_l, seq2.linker_2),
            RightEnd(30),
        )
        .cut(
            sel!(seq2._l_l),
            tr!(seq2._l_l -> seq2._l_l_l, seq2.brc2),
            RightEnd(8),
        )
        .map(
            sel!(seq2.brc2),
            tr!(seq2.brc2 -> seq2.brc2.not_mapped),
            map_file,
            1,
        )
        .cut(
            sel!(seq2._l_l_l),
            tr!(seq2._l_l_l -> seq2._l_l_l_l, seq2.linker_3),
            RightEnd(30),
        )
        .cut(
            sel!(seq2._l_l_l_l),
            tr!(seq2._l_l_l_l -> seq2._l_l_l_l_l, seq2.brc3),
            RightEnd(8),
        )
        .map(
            sel!(seq2.brc3),
            tr!(seq2.brc3 -> seq2.brc3.not_mapped),
            map_file,
            1,
        )
        .cut(
            sel!(seq2._l_l_l_l_l),
            tr!(seq2._l_l_l_l_l -> _, seq2.umi),
            RightEnd(10),
        )
        .dbg(sel!())
        .run()
        .unwrap_or_else(|e| panic!("{e}"));
}
