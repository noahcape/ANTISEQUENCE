use antisequence::*;

fn main() {
    let file1 = "";
    let file2 = "";
    let out1 = "";
    let out2 = "";
    let chunk_size = 256;
    let threads = 6;

    iter_fastq2(file1, file2, chunk_size)
        .unwrap_or_else(|e| panic!("{e}"))
        .match_one(
            sel!(),
            tr!(seq1.* -> seq1.barcode_2, seq1.anchor, seq1._r),
            "CAGAGC",
            HammingSearch(Frac(0.83)),
        )
        .length_in_bounds(
            sel!(seq1.barcode_2),
            tr!(seq1.barcode_2 -> seq1.barcode_2.v_len),
            9..=10,
        )
        .norm(sel!(seq1.barcode_2), label!(seq1.barcode_2), 9..=10)
        .retain(sel!(seq1.barcode_2.v_len, seq1.anchor))
        .cut(
            sel!(seq1._r),
            tr!(seq1._r -> seq1.umi, seq1._r_r),
            LeftEnd(8),
        )
        .retain(sel!(seq1.umi))
        .cut(
            sel!(seq1._r_r),
            tr!(seq1._r_r -> seq1.barcode_2, seq1._r_r_r),
            LeftEnd(10),
        )
        .retain(sel!(seq1.barcode_2))
        .set(sel!(seq2.*), label!(seq2.*), "{seq2.read_2}")
        .set(
            sel!(),
            label!(seq1.*),
            "{seq1.barcode_1}{seq1.umi}{seq1.barcode_2}",
        )
        .collect_fastq2(sel!(), out1, out2)
        .run_with_threads(threads)
}
