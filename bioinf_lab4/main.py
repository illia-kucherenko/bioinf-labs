"""
main.py - Main pipeline for Bioinformatics Lab 4.
Workflow:
  1. Fetch viral genome from NCBI.
  2. Sample 2,000 random reads of 100-120 bp (in-silico shotgun sequencing).
  3. Store reads in an array/list.
  4. Rebuild the DNA sequence using Greedy Overlap Assembly.
  5. Analyze graph topology (De Bruijn) and evaluate against reference.
  6. Diagnose the main algorithmic bottlenecks and limitations.
"""
import json
import os
import sys
from ncbi_fetcher import fetch_ncbi_genome
from sampler import sample_reads, export_reads_fasta
from assembler import GreedyOverlapAssembler, DeBruijnGraphAnalyzer
from evaluator import evaluate_assembly


def plot_diagnostics(stats: dict, assembly_stats: dict, output_dir: str = "data") -> None:
    """Generates visualization plots if matplotlib is available."""
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        print("[Visualizer] matplotlib not installed. Skipping plot generation.")
        return

    os.makedirs(output_dir, exist_ok=True)

    # 1. Coverage Map Plot
    plt.figure(figsize=(10, 4))
    plt.plot(stats["coverage_map"], color="#2b5c8f", lw=1.2)
    plt.axhline(
        y=stats["mean_coverage_depth"],
        color="crimson",
        linestyle="--",
        label=f"Mean Depth ({stats['mean_coverage_depth']}x)",
    )
    plt.title("Genome Sequencing Depth (Coverage Map)", fontsize=13, fontweight="bold")
    plt.xlabel("Genomic Coordinate (bp)", fontsize=11)
    plt.ylabel("Coverage Depth (Number of Reads)", fontsize=11)
    plt.legend()
    plt.tight_layout()
    cov_path = os.path.join(output_dir, "coverage_depth.png")
    plt.savefig(cov_path, dpi=150)
    plt.close()
    print(f"[Visualizer] Coverage plot saved to {cov_path}")

    # 2. Overlap Lengths Distribution
    if assembly_stats["overlap_lengths_distribution"]:
        plt.figure(figsize=(8, 4))
        plt.hist(
            assembly_stats["overlap_lengths_distribution"],
            bins=25,
            color="#3a8a5b",
            edgecolor="black",
        )
        plt.title("Distribution of Overlap Lengths During Assembly", fontsize=13, fontweight="bold")
        plt.xlabel("Overlap Length (bp)", fontsize=11)
        plt.ylabel("Number of Merges", fontsize=11)
        plt.tight_layout()
        overlap_path = os.path.join(output_dir, "overlap_distribution.png")
        plt.savefig(overlap_path, dpi=150)
        plt.close()
        print(f"[Visualizer] Overlap distribution plot saved to {overlap_path}")


def main() -> None:
    print("=" * 70)
    print(" BIOINFORMATICS LAB 4: DE NOVO VIRAL GENOME ASSEMBLY SIMULATION")
    print("=" * 70)

    # Step 1: Sequence Acquisition
    accession = "NC_001422.1"  # Bacteriophage phiX174 (5,386 bp)
    header, reference_genome = fetch_ncbi_genome(accession)
    genome_len = len(reference_genome)

    # Step 2 & 3: In-Silico Shotgun Sampling (2,000 reads, 100-120 bp)
    print("\n" + "-" * 70)
    print("STEP 1: Generating 2,000 Random Reads (100-120 bp)")
    print("-" * 70)
    n_reads = 2000
    min_read_len = 100
    max_read_len = 120

    reads_list, sampling_stats = sample_reads(
        reference_genome,
        n_samples=n_reads,
        min_len=min_read_len,
        max_len=max_read_len,
        seed=42,
    )

    print(f"Total reads stored in list: {len(reads_list):,}")
    print(f"Sample Read 1 ({len(reads_list[0])} bp): {reads_list[0][:60]}...")
    print(f"Total bases sequenced: {sampling_stats['total_sampled_bp']:,} bp")
    print(f"Average Sequencing Depth: {sampling_stats['mean_coverage_depth']}x")
    print(f"Zero-coverage positions: {sampling_stats['zero_coverage_bases']} bp ({sampling_stats['zero_coverage_pct']}%)")

    # Export reads
    export_reads_fasta(reads_list, "data/sampled_reads.fasta")

    # Step 4: De Novo Assembly
    print("\n" + "-" * 70)
    print("STEP 2: Rebuilding DNA Sequence (Greedy Overlap Assembly)")
    print("-" * 70)
    min_overlap = 30
    assembler = GreedyOverlapAssembler(min_overlap=min_overlap)
    contigs, assembly_stats = assembler.assemble(reads_list)

    print(f"Assembly completed in {assembly_stats['assembly_time_seconds']} seconds.")
    print(f"Unique non-contained reads: {assembly_stats['unique_non_contained']}")
    print(f"Total greedy merges performed: {assembly_stats['merges_performed']}")
    print(f"Number of resulting contigs: {assembly_stats['final_contigs_count']}")
    print(f"Longest contig length: {assembly_stats['max_contig_length']:,} bp")

    # Step 5: De Bruijn Graph Analysis
    print("\n" + "-" * 70)
    print("STEP 3: De Bruijn Graph Topological Analysis (k=31)")
    print("-" * 70)
    dbg_analyzer = DeBruijnGraphAnalyzer(k=31)
    dbg_stats = dbg_analyzer.analyze(reads_list)
    print(f"Total (k-1)-mer nodes: {dbg_stats['total_nodes']:,}")
    print(f"Total k-mer edges: {dbg_stats['total_edges']:,}")
    print(f"Branching fork nodes (outdegree > 1): {dbg_stats['branching_fork_nodes']}")
    print(f"Branching join nodes (indegree > 1): {dbg_stats['branching_join_nodes']}")
    print(f"Is graph purely linear (Eulerian path trivial)? {dbg_stats['is_eulerian_linear']}")

    # Step 6: Evaluation against reference
    print("\n" + "-" * 70)
    print("STEP 4: Reference Evaluation & Algorithmic Problem Diagnosis")
    print("-" * 70)
    eval_result = evaluate_assembly(reference_genome, contigs)
    print(f"Reference Length: {eval_result['reference_length']:,} bp")
    print(f"Reconstructed Longest Contig: {eval_result['longest_contig_length']:,} bp")
    print(f"Assembly N50: {eval_result['n50']:,} bp")
    print(f"Genome Coverage Fraction: {eval_result['genome_coverage_fraction']}%")
    print(f"Is Contig an Exact Substring Match? {eval_result['is_exact_match']}")
    if eval_result["is_exact_match"]:
        print(f"Match offset in reference: {eval_result['match_position_in_reference']} bp")

    # Save summary report
    report = {
        "accession": accession,
        "reference_header": header,
        "genome_length": genome_len,
        "sampling_stats": {k: v for k, v in sampling_stats.items() if k != "coverage_map" and k != "metadata"},
        "assembly_stats": {k: v for k, v in assembly_stats.items() if k != "overlap_lengths_distribution"},
        "de_bruijn_analysis": dbg_stats,
        "evaluation": eval_result,
    }

    with open("data/assembly_report.json", "w", encoding="utf-8") as f:
        json.dump(report, f, indent=2)
    print("\n[Report] Saved complete assembly report to data/assembly_report.json")

    # Generate plots
    plot_diagnostics(sampling_stats, assembly_stats)

    # Detailed discussion of the algorithmic bottlenecks
    print("\n" + "=" * 70)
    print(" THE MAIN PROBLEMS WITH THE ALGORITHMIC APPROACH")
    print("=" * 70)
    for idx, prob in enumerate(eval_result["problems_diagnosed"], 1):
        print(f"\n{idx}. [{prob['severity']}] {prob['problem']}")
        print(f"   {prob['impact']}")
    print("=" * 70)


if __name__ == "__main__":
    main()
