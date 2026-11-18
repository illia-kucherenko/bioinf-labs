"""
sampler.py - In-silico shotgun sequencing simulator.
Takes random short-read samples from a target genome sequence.
"""
import random
import math
from typing import TypedDict


class ReadMetadata(TypedDict):
    read_id: int
    start: int
    end: int
    length: int
    sequence: str


def sample_reads(
    genome: str,
    n_samples: int = 2000,
    min_len: int = 100,
    max_len: int = 120,
    seed: int | None = 42,
) -> tuple[list[str], dict]:
    """
    Simulates shotgun sequencing by drawing random substrings from the genome.

    Args:
        genome: The full reference DNA sequence string.
        n_samples: Number of reads to generate (default: 2000).
        min_len: Minimum read length in base pairs (default: 100).
        max_len: Maximum read length in base pairs (default: 120).
        seed: Random seed for reproducibility.

    Returns:
        tuple[list[str], dict]:
            - reads: List of sampled DNA substrings.
            - stats: Detailed coverage and sampling statistics.
    """
    if seed is not None:
        random.seed(seed)

    genome_len = len(genome)
    if genome_len < max_len:
        raise ValueError(f"Genome length ({genome_len}) is shorter than max_len ({max_len}).")

    reads: list[str] = []
    metadata: list[ReadMetadata] = []
    coverage_map = [0] * genome_len

    for i in range(n_samples):
        read_len = random.randint(min_len, max_len)
        start_pos = random.randint(0, genome_len - read_len)
        end_pos = start_pos + read_len
        read_seq = genome[start_pos:end_pos]

        reads.append(read_seq)
        metadata.append({
            "read_id": i,
            "start": start_pos,
            "end": end_pos,
            "length": read_len,
            "sequence": read_seq,
        })

        for pos in range(start_pos, end_pos):
            coverage_map[pos] += 1

    total_sampled_bp = sum(len(r) for r in reads)
    mean_coverage = total_sampled_bp / genome_len
    zero_cov_count = coverage_map.count(0)
    zero_cov_pct = (zero_cov_count / genome_len) * 100

    # Lander-Waterman theoretical zero-coverage expectation
    theoretical_zero_bases = genome_len * math.exp(-mean_coverage)

    stats = {
        "n_samples": n_samples,
        "min_len": min_len,
        "max_len": max_len,
        "total_sampled_bp": total_sampled_bp,
        "genome_length": genome_len,
        "mean_coverage_depth": round(mean_coverage, 2),
        "min_coverage_depth": min(coverage_map),
        "max_coverage_depth": max(coverage_map),
        "zero_coverage_bases": zero_cov_count,
        "zero_coverage_pct": round(zero_cov_pct, 4),
        "lander_waterman_expected_zero": round(theoretical_zero_bases, 2),
        "coverage_map": coverage_map,
        "metadata": metadata,
    }

    return reads, stats


def export_reads_fasta(reads: list[str], output_path: str = "data/reads.fasta") -> None:
    """Exports reads to a standard FASTA file."""
    with open(output_path, "w", encoding="utf-8") as f:
        for idx, read in enumerate(reads):
            f.write(f">read_{idx+1} length={len(read)}\n{read}\n")
    print(f"[Sampler] Exported {len(reads)} reads to {output_path}")


if __name__ == "__main__":
    from ncbi_fetcher import fetch_ncbi_genome
    _, seq = fetch_ncbi_genome("NC_001422.1")
    reads, stats = sample_reads(seq, n_samples=2000, min_len=100, max_len=120)
    print(f"Sampled {len(reads)} reads.")
    print(f"Mean Coverage: {stats['mean_coverage_depth']}x")
    print(f"Zero-coverage bases: {stats['zero_coverage_bases']} ({stats['zero_coverage_pct']}%)")
