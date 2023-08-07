use antisequence::*;

fn main() {
    let file_1 = "./example_data/split_seq_R1.fastq";
    let file_2 = "./example_data/split_seq_R2.fastq";
    let map_file = "./example_data/split_seq_barcodes.txt";

    iter_fastq2(file_1, file_2, 256)
        .unwrap_or_else(|e| panic!("{e}"))
        .cut(sel!(seq2.*), tr!(seq2.* -> seq2.brc1, seq2._r), LeftEnd(8))
        .map(
            sel!(seq2.brc1),
            tr!(seq2.brc1 -> seq2.brc1.not_mapped),
            map_file,
            1,
        )
        .cut(
            sel!(seq2._r),
            tr!(seq2._r -> seq2.linker_2, seq2._r_r),
            LeftEnd(30),
        )
        .cut(
            sel!(seq2._r_r),
            tr!(seq2._r_r -> seq2.brc2, seq2._r_r_r),
            LeftEnd(8),
        )
        .map(
            sel!(seq2.brc2),
            tr!(seq2.brc2 -> seq2.brc2.not_mapped),
            map_file,
            1,
        )
        .cut(
            sel!(seq2._r_r_r),
            tr!(seq2._r_r_r -> seq2.linker_3, seq2._r_r_r_r),
            LeftEnd(30),
        )
        .cut(
            sel!(seq2._r_r_r_r),
            tr!(seq2._r_r_r_r -> seq2.brc3, seq2._r_r_r_r_r),
            LeftEnd(8),
        )
        .map(
            sel!(seq2.brc3),
            tr!(seq2.brc3 -> seq2.brc3.not_mapped),
            map_file,
            1,
        )
        .cut(
            sel!(seq2._r_r_r_r_r),
            tr!(seq2._r_r_r_r_r -> seq2.umi, _),
            LeftEnd(10),
        )
        .dbg(sel!())
        .run()
        .unwrap_or_else(|e| panic!("{e}"));
}
