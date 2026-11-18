"""
evaluator.py - Evaluation and comparative analysis module.
Validates reconstructed contigs against the reference genome and analyzes
the fundamental algorithmic bottlenecks in genome assembly.
"""
from typing import TypedDict


class EvaluationResult(TypedDict):
    reference_length: int
    num_contigs: int
    total_assembled_length: int
    longest_contig_length: int
    n50: int
    genome_coverage_fraction: float
    is_exact_match: bool
    match_position_in_reference: int
    reverse_match: bool
    problems_diagnosed: list[dict]


def calculate_n50(lengths: list[int]) -> int:
    """Calculates the N50 metric for contig lengths."""
    if not lengths:
        return 0
    sorted_lens = sorted(lengths, reverse=True)
    total_length = sum(sorted_lens)
    half_length = total_length / 2.0
    running_sum = 0
    for length in sorted_lens:
        running_sum += length
        if running_sum >= half_length:
            return length
    return sorted_lens[-1]


def evaluate_assembly(reference: str, contigs: list[str]) -> EvaluationResult:
    """
    Evaluates assembled contigs against the reference genome.
    Checks for exact substring matches (including circular wrapping).
    """
    ref_len = len(reference)
    num_contigs = len(contigs)
    lengths = [len(c) for c in contigs]
    total_assembled_len = sum(lengths)
    longest_len = lengths[0] if lengths else 0
    n50 = calculate_n50(lengths)

    # Check match of primary (longest) contig
    primary_contig = contigs[0] if contigs else ""
    # Double reference to account for circular viral genomes (e.g. phiX174)
    ref_doubled = reference + reference

    is_exact = False
    match_pos = -1
    rev_match = False

    if primary_contig:
        pos = ref_doubled.find(primary_contig)
        if pos != -1:
            is_exact = True
            match_pos = pos % ref_len
        else:
            # Check reverse complement
            trans = str.maketrans("ACGTN", "TGCAN")
            rc_primary = primary_contig.translate(trans)[::-1]
            pos_rc = ref_doubled.find(rc_primary)
            if pos_rc != -1:
                is_exact = True
                match_pos = pos_rc % ref_len
                rev_match = True

    # Genome fraction covered
    cov_mask = [False] * ref_len
    for c in contigs:
        pos = ref_doubled.find(c)
        if pos != -1:
            for p in range(pos, pos + len(c)):
                cov_mask[p % ref_len] = True

    covered_bases = sum(1 for b in cov_mask if b)
    coverage_fraction = (covered_bases / ref_len) * 100

    # Algorithmic Problems Diagnosis
    problems = [
        {
            "problem": "Repetitive Sequences (The Assembly Bottleneck)",
            "impact": (
                "When identical repeats are longer than the overlap cutoff, the algorithm cannot "
                "determine the correct linear order of reads flanking the repeat. "
                "Greedy assemblers will either collapse the repeated copies into a single instance "
                "or create false chimeric joins."
            ),
            "severity": "CRITICAL"
        },
        {
            "problem": "Coverage Gaps (Lander-Waterman Dropout)",
            "impact": (
                "Random sampling leaves regions with zero reads (or insufficient overlap < K). "
                f"Even with 2,000 reads, any gap breaks the assembly into {num_contigs} separate "
                "disconnected contigs rather than a single chromosome."
            ),
            "severity": "HIGH" if num_contigs > 1 else "LOW"
        },
        {
            "problem": "Computational Complexity (NP-Hard Shortest Superstring)",
            "impact": (
                "Reconstructing the minimal sequence that contains all reads is the NP-hard "
                "Shortest Superstring Problem. Finding all pairwise overlaps naively requires O(N^2 * L) "
                "comparisons. Greedy heuristics get trapped in local optima and cannot backtrack."
            ),
            "severity": "HIGH"
        },
        {
            "problem": "Circularity and Boundary Overhangs",
            "impact": (
                "Viral genomes (like PhiX174) are circular molecules. Linear shotgun reads that span the "
                "origin (0 / 5,386 bp) cannot be linearly assembled unless the assembler handles circular "
                "boundary conditions, leaving terminal overhangs."
            ),
            "severity": "MEDIUM"
        }
    ]

    return {
        "reference_length": ref_len,
        "num_contigs": num_contigs,
        "total_assembled_length": total_assembled_len,
        "longest_contig_length": longest_len,
        "n50": n50,
        "genome_coverage_fraction": round(coverage_fraction, 2),
        "is_exact_match": is_exact,
        "match_position_in_reference": match_pos,
        "reverse_match": rev_match,
        "problems_diagnosed": problems,
    }
