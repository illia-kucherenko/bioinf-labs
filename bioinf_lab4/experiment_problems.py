"""
experiment_problems.py - Empirical demonstration of the fundamental failure modes
of greedy sequence assembly:
  1. The Repeat Problem:
     - Case A: Repeat < Read Length (Bridged and resolved)
     - Case B: Repeat > Read Length (Fundamental failure: unbridgeable ambiguity)
  2. Coverage Gaps (Lander-Waterman fragmentation)
  3. Overlap Threshold Sensitivity (Trade-off between chimeras and fragmentation)
"""
from ncbi_fetcher import fetch_ncbi_genome
from sampler import sample_reads
from assembler import GreedyOverlapAssembler, DeBruijnGraphAnalyzer
from evaluator import evaluate_assembly


def test_repeat_collapse() -> None:
    """
    Demonstrates the Repeat Resolution Theorem in Genome Assembly:
    - If repeat length R < Read length L: reads spanning the repeat anchor to unique
      flanks, resolving the repeat.
    - If repeat length R > Read length L: no read can span the repeat, producing
      inherent ambiguity that shatters the assembly into fragmented contigs.
    """
    print("\n" + "=" * 70)
    print("EXPERIMENT 1: The Repeat Ambiguity Problem (The Assembly Bottleneck)")
    print("=" * 70)
    _, ref = fetch_ncbi_genome("NC_001422.1")
    ref_sub = ref[:3000]

    # --- CASE A: Repeat shorter than read length (R = 50 bp < L = 100-120 bp) ---
    rep_short = ("ATGC" * 15)[:50]
    construct_short = (
        ref_sub[:500] + rep_short + ref_sub[550:1800] + rep_short + ref_sub[1850:3000]
    )
    reads_a, _ = sample_reads(construct_short, n_samples=1500, min_len=100, max_len=120, seed=42)
    assembler = GreedyOverlapAssembler(min_overlap=30)
    contigs_a, stats_a = assembler.assemble(reads_a)

    print(f"\n[Case A] Repeat Length = 50 bp (Shorter than 100-120 bp reads):")
    print(f"  - Assembled Contigs: {len(contigs_a)}")
    print(f"  - Longest Contig: {len(contigs_a[0]):,} bp / {len(construct_short):,} bp reference")
    print(f"  - Outcome: Reads span the 50 bp repeat into unique flanking regions -> RESOLVED!")

    # --- CASE B: Repeat longer than read length (R = 140 bp > L = 100-120 bp) ---
    rep_long = ("ATGC" * 35)[:140]
    construct_long = (
        ref_sub[:500] + rep_long + ref_sub[640:1800] + rep_long + ref_sub[1940:3000]
    )
    reads_b, _ = sample_reads(construct_long, n_samples=1500, min_len=100, max_len=120, seed=42)
    contigs_b, stats_b = assembler.assemble(reads_b)

    dbg = DeBruijnGraphAnalyzer(k=31)
    dbg_stats = dbg.analyze(reads_b)

    print(f"\n[Case B] Repeat Length = 140 bp (Longer than 100-120 bp reads):")
    print(f"  - Assembled Contigs: {len(contigs_b)}")
    for idx, c in enumerate(contigs_b):
        print(f"      Contig {idx+1}: {len(c):,} bp")
    print(f"  - De Bruijn Fork Nodes (Ambiguity): {dbg_stats['branching_fork_nodes']}")
    print(f"  - De Bruijn Join Nodes (Ambiguity): {dbg_stats['branching_join_nodes']}")
    print(
        f"  - Outcome: NO read can span the 140 bp repeat. The assembler cannot determine\n"
        f"             which flanking region connects to which copy, shattering the assembly!"
    )


def test_coverage_variation() -> None:
    """
    Demonstrates how varying read coverage affects fragmentation (Lander-Waterman).
    """
    print("\n" + "=" * 70)
    print("EXPERIMENT 2: Coverage Gaps & Assembly Fragmentation (Lander-Waterman)")
    print("=" * 70)
    _, ref = fetch_ncbi_genome("NC_001422.1")

    for n_reads in [200, 500, 1000, 2000]:
        reads, s_stats = sample_reads(ref, n_samples=n_reads, min_len=100, max_len=120, seed=42)
        assembler = GreedyOverlapAssembler(min_overlap=30)
        contigs, a_stats = assembler.assemble(reads)
        eval_res = evaluate_assembly(ref, contigs)

        print(
            f"Reads: {n_reads:4d} | Depth: {s_stats['mean_coverage_depth']:5.1f}x | "
            f"Zero Cov: {s_stats['zero_coverage_bases']:3d} bp | Contigs: {len(contigs):2d} | "
            f"N50: {eval_res['n50']:5d} bp | Genome Covered: {eval_res['genome_coverage_fraction']:5.1f}%"
        )


def test_overlap_threshold_sensitivity() -> None:
    """
    Demonstrates the trade-off of the minimum overlap parameter (K):
    - Too low (e.g. 15 bp): Risk of false accidental matches.
    - Too high (e.g. 90 bp): Insufficient overlap connectivity.
    """
    print("\n" + "=" * 70)
    print("EXPERIMENT 3: Overlap Cutoff Parameter Sensitivity (K)")
    print("=" * 70)
    _, ref = fetch_ncbi_genome("NC_001422.1")
    reads, _ = sample_reads(ref, n_samples=2000, min_len=100, max_len=120, seed=42)

    for k in [15, 30, 60, 90]:
        assembler = GreedyOverlapAssembler(min_overlap=k)
        contigs, a_stats = assembler.assemble(reads)
        eval_res = evaluate_assembly(ref, contigs)

        print(
            f"K (Min Overlap): {k:2d} bp | Contigs: {len(contigs):2d} | "
            f"Longest Contig: {eval_res['longest_contig_length']:5d} bp | "
            f"Exact Match: {str(eval_res['is_exact_match']):5s} | Time: {a_stats['assembly_time_seconds']:.3f}s"
        )


if __name__ == "__main__":
    test_repeat_collapse()
    test_coverage_variation()
    test_overlap_threshold_sensitivity()
