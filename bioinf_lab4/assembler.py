"""
assembler.py - De novo sequence assembly using Greedy Overlap (Shortest Superstring)
and De Bruijn Graph topological analysis.
"""
import heapq
import time
from collections import defaultdict
from typing import TypedDict


class AssemblyStats(TypedDict):
    initial_reads: int
    unique_non_contained: int
    merges_performed: int
    final_contigs_count: int
    max_contig_length: int
    min_contig_length: int
    contig_lengths: list[int]
    assembly_time_seconds: float
    min_overlap_threshold: int
    ambiguous_tie_merges: int
    overlap_lengths_distribution: list[int]


def remove_contained_reads(reads: list[str]) -> list[str]:
    """
    Filters out exact duplicate reads and reads that are substrings
    of longer reads.
    """
    unique_reads = sorted(list(set(reads)), key=len, reverse=True)
    filtered: list[str] = []
    for r in unique_reads:
        # Check if r is already a substring of any keep read
        is_contained = False
        for parent in filtered:
            if r in parent:
                is_contained = True
                break
        if not is_contained:
            filtered.append(r)
    return filtered


class GreedyOverlapAssembler:
    """
    Greedy Shortest Superstring Assembler using an inverted prefix index
    and max-heap for fast O(N log N) overlap resolution.
    """

    def __init__(self, min_overlap: int = 30):
        self.min_overlap = min_overlap

    @staticmethod
    def _find_overlap(s1: str, s2: str, min_len: int) -> int:
        """Finds the maximum suffix of s1 that is a prefix of s2 (length >= min_len)."""
        start = 0
        while True:
            pos = s1.find(s2[:min_len], start)
            if pos == -1:
                return 0
            if s2.startswith(s1[pos:]):
                return len(s1) - pos
            start = pos + 1

    def assemble(self, reads: list[str]) -> tuple[list[str], AssemblyStats]:
        """
        Reconstructs sequence contigs from reads greedily merging highest overlapping pairs.
        """
        start_time = time.time()
        initial_count = len(reads)

        # 1. Deduplicate & remove contained reads
        filtered_reads = remove_contained_reads(reads)
        unique_count = len(filtered_reads)

        # Active contigs indexed by integer ID
        contigs: dict[int, str] = {i: seq for i, seq in enumerate(filtered_reads)}
        next_id = unique_count

        # 2. Build prefix index for fast candidate retrieval
        k = self.min_overlap
        prefix_to_ids: dict[str, set[int]] = defaultdict(set)
        for cid, seq in contigs.items():
            if len(seq) >= k:
                prefix_to_ids[seq[:k]].add(cid)

        # 3. Find candidate overlaps & populate max-heap (-overlap_len, id1, id2)
        # We store negative overlap for min-heap implementation of max-heap
        heap: list[tuple[int, int, int]] = []
        ambiguous_ties = 0
        overlap_dist: list[int] = []

        def get_overlaps_for(cid: int, seq: str):
            """Finds overlaps where seq is the left (suffix) partner."""
            results = []
            if len(seq) < k:
                return results
            # Suffixes must match prefix of another contig
            for pos in range(len(seq) - k + 1):
                kmer = seq[pos:pos + k]
                if kmer in prefix_to_ids:
                    suffix = seq[pos:]
                    suffix_len = len(suffix)
                    for target_id in prefix_to_ids[kmer]:
                        if target_id != cid and contigs[target_id].startswith(suffix):
                            results.append((-suffix_len, cid, target_id))
            return results

        # Initial all-pairs overlap scan
        for cid, seq in contigs.items():
            for entry in get_overlaps_for(cid, seq):
                heapq.heappush(heap, entry)

        merges = 0

        # 4. Iterative Greedy Merging
        while heap:
            neg_olen, id1, id2 = heapq.heappop(heap)
            olen = -neg_olen

            if olen < self.min_overlap:
                break

            # Check if both contigs are still active
            if id1 not in contigs or id2 not in contigs:
                continue

            seq1 = contigs[id1]
            seq2 = contigs[id2]

            # Verify overlap is still valid
            if not seq1.endswith(seq2[:olen]):
                continue

            # Check for ambiguity: are there other merges with identical score?
            if heap and -heap[0][0] == olen and (heap[0][1] == id1 or heap[0][2] == id2):
                ambiguous_ties += 1

            # Merge: seq1 + remaining part of seq2
            merged_seq = seq1 + seq2[olen:]
            del contigs[id1]
            del contigs[id2]

            # Update prefix index
            prefix_to_ids[seq1[:k]].discard(id1)
            prefix_to_ids[seq2[:k]].discard(id2)

            new_cid = next_id
            next_id += 1
            contigs[new_cid] = merged_seq

            if len(merged_seq) >= k:
                prefix_to_ids[merged_seq[:k]].add(new_cid)

            # Compute new overlaps for the merged contig
            # (as suffix partner)
            for entry in get_overlaps_for(new_cid, merged_seq):
                heapq.heappush(heap, entry)

            # (as prefix partner for remaining active contigs)
            new_prefix = merged_seq[:k]
            for other_id, other_seq in contigs.items():
                if other_id == new_cid or len(other_seq) < k:
                    continue
                pos = other_seq.find(new_prefix)
                while pos != -1:
                    suffix = other_seq[pos:]
                    if merged_seq.startswith(suffix) and len(suffix) >= self.min_overlap:
                        heapq.heappush(heap, (-len(suffix), other_id, new_cid))
                    pos = other_seq.find(new_prefix, pos + 1)

            merges += 1
            overlap_dist.append(olen)

        final_contigs = sorted(list(contigs.values()), key=len, reverse=True)
        elapsed = time.time() - start_time
        contig_lengths = [len(c) for c in final_contigs]

        stats: AssemblyStats = {
            "initial_reads": initial_count,
            "unique_non_contained": unique_count,
            "merges_performed": merges,
            "final_contigs_count": len(final_contigs),
            "max_contig_length": contig_lengths[0] if contig_lengths else 0,
            "min_contig_length": contig_lengths[-1] if contig_lengths else 0,
            "contig_lengths": contig_lengths,
            "assembly_time_seconds": round(elapsed, 3),
            "min_overlap_threshold": self.min_overlap,
            "ambiguous_tie_merges": ambiguous_ties,
            "overlap_lengths_distribution": overlap_dist,
        }

        return final_contigs, stats


class DeBruijnGraphAnalyzer:
    """
    Decomposes reads into (k)-mers to construct a De Bruijn Graph and analyze
    branching nodes (indegree > 1 or outdegree > 1), demonstrating why repeats
    shatter sequence assemblies.
    """

    def __init__(self, k: int = 31):
        self.k = k

    def analyze(self, reads: list[str]) -> dict:
        """
        Builds the De Bruijn graph adjacency from reads and detects repeat bifurcations.
        """
        graph: dict[str, set[str]] = defaultdict(set)
        in_degrees: dict[str, int] = defaultdict(int)
        out_degrees: dict[str, int] = defaultdict(int)

        for read in reads:
            for i in range(len(read) - self.k + 1):
                kmer = read[i:i + self.k]
                u = kmer[:-1]  # Prefix (k-1)-mer
                v = kmer[1:]   # Suffix (k-1)-mer
                if v not in graph[u]:
                    graph[u].add(v)
                    out_degrees[u] += 1
                    in_degrees[v] += 1

        all_nodes = set(graph.keys()) | set(in_degrees.keys())
        fork_nodes = [n for n in all_nodes if out_degrees[n] > 1]
        join_nodes = [n for n in all_nodes if in_degrees[n] > 1]
        complex_nodes = [n for n in all_nodes if out_degrees[n] > 1 and in_degrees[n] > 1]

        return {
            "k": self.k,
            "total_nodes": len(all_nodes),
            "total_edges": sum(len(neighbors) for neighbors in graph.values()),
            "branching_fork_nodes": len(fork_nodes),
            "branching_join_nodes": len(join_nodes),
            "complex_nodes": len(complex_nodes),
            "is_eulerian_linear": len(fork_nodes) == 0 and len(join_nodes) == 0,
        }
