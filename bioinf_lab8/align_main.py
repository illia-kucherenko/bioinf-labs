"""Runner script for Task 2: Pairwise DNA Sequence Alignment using Matrices."""

import os
from sequence_alignment import needleman_wunsch, format_matrix_text, plot_alignment_matrix


def run_alignment_demo():
    output_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "output")
    os.makedirs(output_dir, exist_ok=True)

    # Example 1: Classic benchmark pair
    seq1 = "GATTACA"
    seq2 = "GCATGCA"
    match = 1
    mismatch = -1
    gap = -1

    print("Task 2: Pairwise DNA Sequence Alignment using Matrices")
    print(f"Sequence 1: {seq1}")
    print(f"Sequence 2: {seq2}")
    print(f"Parameters: Match = {match:+d}, Mismatch = {mismatch:+d}, Gap = {gap:+d}\n")

    res = needleman_wunsch(seq1, seq2, match=match, mismatch=mismatch, gap=gap)

    print("Dynamic Programming Score Matrix (* indicates optimal traceback path):")
    matrix_str = format_matrix_text(res["score_matrix"], seq1, seq2, res["path"])
    print(matrix_str)

    print("\nOptimal Global Alignment:")
    print(f"  Seq 1:   {res['aligned_seq1']}")
    print(f"           {res['alignment_marks']}")
    print(f"  Seq 2:   {res['aligned_seq2']}")

    stats = res["stats"]
    print(f"\nAlignment Score: {res['score']}")
    print(f"Matches:    {stats['matches']}")
    print(f"Mismatches: {stats['mismatches']}")
    print(f"Gaps:       {stats['gaps']}")
    print(f"Identity:   {stats['identity']:.1f}%")

    plot_path = os.path.join(output_dir, "04_alignment_matrix.png")
    plot_alignment_matrix(
        res["score_matrix"],
        seq1,
        seq2,
        path=res["path"],
        title=f"Needleman-Wunsch Alignment ({seq1} vs {seq2})",
        save_path=plot_path
    )

    # Example 2: Short motif pair with gaps
    seqA = "AGCACG"
    seqB = "ACAC"
    res2 = needleman_wunsch(seqA, seqB, match=2, mismatch=-1, gap=-2)
    print(f"\nExample 2: {seqA} vs {seqB} (Match=+2, Mismatch=-1, Gap=-2)")
    print(f"  Seq A:   {res2['aligned_seq1']}")
    print(f"           {res2['alignment_marks']}")
    print(f"  Seq B:   {res2['aligned_seq2']}")
    print(f"  Score:   {res2['score']}")

    plot_path_2 = os.path.join(output_dir, "05_alignment_matrix_example2.png")
    plot_alignment_matrix(
        res2["score_matrix"],
        seqA,
        seqB,
        path=res2["path"],
        title=f"Alignment Matrix ({seqA} vs {seqB})",
        save_path=plot_path_2
    )


if __name__ == "__main__":
    run_alignment_demo()
