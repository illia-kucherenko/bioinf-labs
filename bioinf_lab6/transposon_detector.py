"""
Transposon Detector Module
Bioinformatics Lab 6 - Detection of Transposable Elements in DNA Sequences
Implements both Homology-based (Seed-and-Extend + Local Alignment + Fragment Chaining)
and Structural (De Novo TIR & TSD) detection with configurable in-between spacing.
"""

from dataclasses import dataclass, field
from typing import List, Dict, Tuple, Optional, Set
from transposon_simulator import reverse_complement, TransposableElement


@dataclass
class DetectionResult:
    """Represents a detected transposable element or fragment."""
    element_name: str
    start: int          # 0-based inclusive
    end: int            # 0-based inclusive
    method: str         # "Homology (Local Alignment)", "Structural (De Novo TIR)", "Chained Nested Span", "New Candidate TE"
    score: float        # Alignment score or repeat quality score
    identity: float     # Percentage identity (0-100%)
    status: str         # "INTACT", "FRAGMENT_5P", "FRAGMENT_3P", "INTERNAL_FRAGMENT", "CHAINED_INTERRUPTED", "NEW_CANDIDATE_TIR", "STRUCTURAL_NESTED_PARENT"
    query_start: int = 0
    query_end: int = 0
    query_length: int = 0
    in_between_len: int = 0
    left_tir: str = ""
    right_tir: str = ""
    tsd: str = ""
    interrupted_by: List[str] = field(default_factory=list)
    nested_in: Optional[str] = None
    details: str = ""

    @property
    def length(self) -> int:
        return self.end - self.start + 1

    @property
    def start_1based(self) -> int:
        return self.start + 1

    @property
    def end_1based(self) -> int:
        return self.end + 1

    @property
    def coverage(self) -> float:
        """Percentage of the query element covered by this alignment."""
        if self.query_length <= 0:
            return 100.0
        return min(100.0, (self.query_end - self.query_start + 1) / self.query_length * 100.0)


def find_maximal_exact_matches(
    genome: str,
    query: str,
    min_len: int = 12
) -> List[Tuple[int, int, int, int, int]]:
    """
    Find maximal exact matches between genome and query sequence.
    Returns: list of (g_start, g_end, q_start, q_end, length)
    """
    matches = []
    k = min_len
    seen: Set[Tuple[int, int, int, int]] = set()

    for q_start in range(len(query) - k + 1):
        kmer = query[q_start : q_start + k]
        idx = 0
        while True:
            pos = genome.find(kmer, idx)
            if pos == -1:
                break

            # Extend right
            ext_len = k
            while (q_start + ext_len < len(query) and
                   pos + ext_len < len(genome) and
                   query[q_start + ext_len] == genome[pos + ext_len]):
                ext_len += 1

            # Extend left
            left_ext = 0
            while (q_start - left_ext > 0 and
                   pos - left_ext > 0 and
                   query[q_start - left_ext - 1] == genome[pos - left_ext - 1]):
                left_ext += 1

            g_start = pos - left_ext
            g_end = pos + ext_len - 1
            q_s = q_start - left_ext
            q_e = q_start + ext_len - 1
            match_len = g_end - g_start + 1

            key = (g_start, g_end, q_s, q_e)
            if key not in seen:
                seen.add(key)
                matches.append((g_start, g_end, q_s, q_e, match_len))

            idx = pos + 1

    # Filter non-maximal matches
    matches.sort(key=lambda x: -x[4])
    maximal = []
    for m in matches:
        g1, g2, q1, q2, l = m
        if any(g1 >= ex[0] and g2 <= ex[1] and q1 >= ex[2] and q2 <= ex[3] for ex in maximal):
            continue
        maximal.append(m)

    maximal.sort(key=lambda x: x[0])
    return maximal


def smith_waterman(
    seq1: str,
    seq2: str,
    match_score: int = 2,
    mismatch_penalty: int = -2,
    gap_open: int = -6,
    gap_extend: int = -2
) -> List[Tuple[int, int, int, int, float, float, str, str]]:
    """
    Smith-Waterman local alignment for detecting mutated or gapped regions.
    seq1: Genome (target)
    seq2: Transposon (query)
    """
    n, m = len(seq1), len(seq2)
    if n == 0 or m == 0:
        return []

    H = [[0] * (m + 1) for _ in range(n + 1)]
    TB = [[0] * (m + 1) for _ in range(n + 1)]

    for i in range(1, n + 1):
        c1 = seq1[i - 1]
        for j in range(1, m + 1):
            c2 = seq2[j - 1]
            diag = H[i - 1][j - 1] + (match_score if c1 == c2 else mismatch_penalty)
            up = H[i - 1][j] + gap_extend
            left = H[i][j - 1] + gap_extend

            val = max(0, diag, up, left)
            H[i][j] = val

            if val == 0:
                TB[i][j] = 0
            elif val == diag:
                TB[i][j] = 1
            elif val == up:
                TB[i][j] = 2
            else:
                TB[i][j] = 3

    alignments = []
    used_cells = set()

    candidates = []
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            if H[i][j] >= 16:
                candidates.append((H[i][j], i, j))
    candidates.sort(reverse=True, key=lambda x: x[0])

    for score, i_curr, j_curr in candidates:
        if (i_curr, j_curr) in used_cells:
            continue

        i, j = i_curr, j_curr
        align1, align2, path = [], [], []
        matches = 0
        total_len = 0

        while i > 0 and j > 0 and H[i][j] > 0:
            path.append((i, j))
            direction = TB[i][j]
            if direction == 1:
                c1 = seq1[i - 1]
                c2 = seq2[j - 1]
                align1.append(c1)
                align2.append(c2)
                if c1 == c2:
                    matches += 1
                total_len += 1
                i -= 1
                j -= 1
            elif direction == 2:
                align1.append(seq1[i - 1])
                align2.append('-')
                total_len += 1
                i -= 1
            elif direction == 3:
                align1.append('-')
                align2.append(seq2[j - 1])
                total_len += 1
                j -= 1
            else:
                break

        if total_len < 10:
            continue

        path_set = set(path)
        if len(path_set.intersection(used_cells)) > len(path_set) * 0.25:
            continue

        used_cells.update(path_set)

        align1_str = "".join(reversed(align1))
        align2_str = "".join(reversed(align2))
        s1_start = i
        s1_end = i_curr - 1
        s2_start = j
        s2_end = j_curr - 1
        ident = (matches / total_len) * 100.0 if total_len > 0 else 0.0

        alignments.append((s1_start, s1_end, s2_start, s2_end, float(score), ident, align1_str, align2_str))

    return alignments


class HomologyDetector:
    """
    Detects transposable elements via sequence homology to a TE library,
    capable of identifying intact elements, split fragments, and chained nested loci.
    """

    def __init__(self, min_fragment_len: int = 12, min_identity: float = 80.0):
        self.min_fragment_len = min_fragment_len
        self.min_identity = min_identity

    def detect(
        self,
        genome_seq: str,
        library: Dict[str, TransposableElement]
    ) -> List[DetectionResult]:
        """
        Scan genome against each element in library and detect full/fragmented occurrences.
        """
        raw_detections: List[DetectionResult] = []

        for name, te in library.items():
            exact_matches = find_maximal_exact_matches(genome_seq, te.sequence, min_len=self.min_fragment_len)
            covered_regions = set()

            for g_start, g_end, q_start, q_end, match_len in exact_matches:
                covered_regions.add((g_start, g_end))
                cov = (q_end - q_start + 1) / te.length * 100.0
                is_5p = (q_start <= 1)
                is_3p = (q_end >= te.length - 2)

                if is_5p and is_3p and cov >= 85.0:
                    status = "INTACT"
                elif is_5p and not is_3p:
                    status = "FRAGMENT_5P"
                elif is_3p and not is_5p:
                    status = "FRAGMENT_3P"
                else:
                    status = "INTERNAL_FRAGMENT"

                raw_detections.append(DetectionResult(
                    element_name=name,
                    start=g_start,
                    end=g_end,
                    method="Homology (Local Alignment)",
                    score=float(match_len * 2),
                    identity=100.0,
                    status=status,
                    query_start=q_start,
                    query_end=q_end,
                    query_length=te.length,
                    details=f"Exact match [{q_start}:{q_end}], length={match_len}bp, Cov={cov:.1f}%"
                ))

            # Run Smith-Waterman to pick up any mutated occurrences
            sw_alignments = smith_waterman(genome_seq, te.sequence)
            for s1_start, s1_end, s2_start, s2_end, score, ident, a1, a2 in sw_alignments:
                frag_len = s1_end - s1_start + 1
                if frag_len < self.min_fragment_len or ident < self.min_identity:
                    continue

                if any(abs(s1_start - g1) <= 2 and abs(s1_end - g2) <= 2 for g1, g2 in covered_regions):
                    continue

                cov = (s2_end - s2_start + 1) / te.length * 100.0
                is_5p = (s2_start <= 2)
                is_3p = (s2_end >= te.length - 3)

                if is_5p and is_3p and cov >= 80.0:
                    status = "INTACT"
                elif is_5p and not is_3p:
                    status = "FRAGMENT_5P"
                elif is_3p and not is_5p:
                    status = "FRAGMENT_3P"
                else:
                    status = "INTERNAL_FRAGMENT"

                raw_detections.append(DetectionResult(
                    element_name=name,
                    start=s1_start,
                    end=s1_end,
                    method="Homology (Local Alignment)",
                    score=score,
                    identity=ident,
                    status=status,
                    query_start=s2_start,
                    query_end=s2_end,
                    query_length=te.length,
                    details=f"Query[{s2_start}:{s2_end}], Match Len={frag_len}bp, Ident={ident:.1f}%, Cov={cov:.1f}%"
                ))

        raw_detections.sort(key=lambda d: d.start)
        chained_detections = self._chain_interrupted_elements(raw_detections, library)
        return chained_detections

    def _chain_interrupted_elements(
        self,
        detections: List[DetectionResult],
        library: Dict[str, TransposableElement]
    ) -> List[DetectionResult]:
        """
        Detect nested insertion interruptions:
        e.g., TE1_5P ... [TE2] ... TE1_3P
        Returns atomic fragments plus chained compound locus spans.
        """
        final_list = list(detections)

        by_element: Dict[str, List[DetectionResult]] = {}
        for d in detections:
            by_element.setdefault(d.element_name, []).append(d)

        for el_name, frags in by_element.items():
            if len(frags) >= 2:
                frags.sort(key=lambda x: x.start)
                frag_5p = next((f for f in frags if f.status == "FRAGMENT_5P"), None)
                frag_3p = next((f for f in frags if f.status == "FRAGMENT_3P"), None)

                if frag_5p and frag_3p and frag_5p.end < frag_3p.start:
                    intervening_start = frag_5p.end + 1
                    intervening_end = frag_3p.start - 1
                    intervening_len = intervening_end - intervening_start + 1

                    interrupters = []
                    for other in detections:
                        if other.element_name != el_name:
                            if other.start >= frag_5p.end - 2 and other.end <= frag_3p.start + 2:
                                interrupters.append(other.element_name)

                    interrupters = list(dict.fromkeys(interrupters))

                    frag_5p.interrupted_by = interrupters
                    frag_3p.interrupted_by = interrupters

                    combined_score = frag_5p.score + frag_3p.score
                    chained = DetectionResult(
                        element_name=f"{el_name} [Chained Interrupted Locus]",
                        start=frag_5p.start,
                        end=frag_3p.end,
                        method="Chained Nested Span",
                        score=combined_score,
                        identity=(frag_5p.identity + frag_3p.identity) / 2.0,
                        status="CHAINED_INTERRUPTED",
                        query_start=frag_5p.query_start,
                        query_end=frag_3p.query_end,
                        query_length=library[el_name].length,
                        in_between_len=intervening_len,
                        interrupted_by=interrupters,
                        details=(
                            f"Split into 2 fragments across [{intervening_start}:{intervening_end}] "
                            f"({intervening_len} bp in-between). Interrupted by: {', '.join(interrupters) if interrupters else 'Nested TE'}."
                        )
                    )
                    final_list.append(chained)

        final_list.sort(key=lambda d: (d.start, -d.length))
        return final_list


class StructuralDetector:
    """
    De Novo Transposon Detector based on structural features:
    Finds Terminal Inverted Repeats (TIRs) and Target Site Duplications (TSDs)
    without prior knowledge of TE sequences.
    
    Supports strict or minimized spacing in between the two TIRs (default: min_between=10 bp).
    """

    def __init__(
        self,
        min_tir_length: int = 8,
        max_tir_length: int = 16,
        min_between: int = 10,       # Minimum number of bases in between left and right TIRs
        max_between: int = 400,      # Maximum number of bases in between left and right TIRs
        max_mismatches: int = 0
    ):
        self.min_tir_length = min_tir_length
        self.max_tir_length = max_tir_length
        self.min_between = max(1, min_between)
        self.max_between = max_between
        self.max_mismatches = max_mismatches

    def _hamming_distance(self, s1: str, s2: str) -> int:
        return sum(1 for a, b in zip(s1, s2) if a != b)

    def detect(
        self,
        genome_seq: str,
        min_between: Optional[int] = None,
        min_tir_length: Optional[int] = None
    ) -> List[DetectionResult]:
        """
        Scan DNA sequence for pairs of inverted repeats forming putative new transposons.
        The internal region in between the TIRs is guaranteed to be >= min_between bases.
        """
        if min_between is not None:
            min_between = max(1, min_between)
        else:
            min_between = self.min_between

        if min_tir_length is not None:
            tir_min = max(4, min_tir_length)
        else:
            tir_min = self.min_tir_length

        n = len(genome_seq)
        results: List[DetectionResult] = []

        # Scan across TIR lengths
        for tir_len in range(self.max_tir_length, tir_min - 1, -1):
            min_element_len = 2 * tir_len + min_between
            if n < min_element_len:
                continue

            for i in range(0, n - min_element_len + 1):
                left_tir = genome_seq[i : i + tir_len]
                expected_right_tir = reverse_complement(left_tir)

                # Ensure distance in between TIRs is at least min_between bases:
                # Left TIR ends at i + tir_len - 1.
                # Right TIR starts at j.
                # In-between bases = j - (i + tir_len) >= min_between  ==> j >= i + tir_len + min_between
                min_j = i + tir_len + min_between
                max_j = min(n - tir_len, i + tir_len + self.max_between)

                for j in range(min_j, max_j + 1):
                    right_tir = genome_seq[j : j + tir_len]
                    dist = self._hamming_distance(right_tir, expected_right_tir)

                    if dist <= self.max_mismatches:
                        te_start = i
                        te_end = j + tir_len - 1
                        te_len = te_end - te_start + 1
                        actual_between_len = j - (i + tir_len)

                        # Check for Target Site Duplication (TSD) flanking the TIRs (2, 3, 4 bp)
                        tsd_found = False
                        tsd_seq = ""
                        for tsd_len in [4, 3, 2]:
                            if te_start >= tsd_len and te_end + tsd_len < n:
                                tsd_left = genome_seq[te_start - tsd_len : te_start]
                                tsd_right = genome_seq[te_end + 1 : te_end + 1 + tsd_len]
                                if tsd_left == tsd_right:
                                    tsd_found = True
                                    tsd_seq = tsd_left
                                    break

                        score = (tir_len * 10) - (dist * 15) + (25 if tsd_found else 0)
                        ident = ((tir_len - dist) / tir_len) * 100.0

                        details = f"TIR={left_tir} (len {tir_len}bp), In-between={actual_between_len}bp"
                        if tsd_found:
                            details += f", TSD={tsd_seq} ({len(tsd_seq)}bp)"

                        results.append(DetectionResult(
                            element_name=f"New_Putative_TE_{te_start}_{te_end}",
                            start=te_start,
                            end=te_end,
                            method="Structural (De Novo TIR)",
                            score=float(score),
                            identity=ident,
                            status="NEW_CANDIDATE_TIR",
                            query_start=0,
                            query_end=te_len - 1,
                            query_length=te_len,
                            in_between_len=actual_between_len,
                            left_tir=left_tir,
                            right_tir=right_tir,
                            tsd=tsd_seq,
                            details=details
                        ))

        # Filter redundant overlapping candidates, keeping highest scoring non-redundant representatives
        filtered: List[DetectionResult] = []
        results.sort(key=lambda r: (-r.score, -(r.end - r.start), r.start))

        for cand in results:
            overlap = False
            for ex in filtered:
                # Check if almost identical interval (IoU > 0.65 or boundaries within 4bp)
                ov = max(0, min(cand.end, ex.end) - max(cand.start, ex.start) + 1)
                union = cand.length + ex.length - ov
                iou = ov / union if union > 0 else 0.0

                if iou > 0.65 or (abs(cand.start - ex.start) <= 4 and abs(cand.end - ex.end) <= 4):
                    overlap = True
                    break
            if not overlap:
                filtered.append(cand)

        # Detect structural nesting (an inner TIR pair located strictly inside an outer TIR pair)
        for cand in filtered:
            inner_elements = []
            for other in filtered:
                if other != cand and other.start > cand.start + 5 and other.end < cand.end - 5:
                    inner_elements.append(other.element_name)
            if inner_elements:
                cand.status = "STRUCTURAL_NESTED_PARENT"
                cand.interrupted_by = inner_elements
                cand.details += f" [Encloses nested candidate(s): {', '.join(inner_elements)}]"

        filtered.sort(key=lambda r: r.start)
        return filtered


class TransposonDetector:
    """
    Integrated Transposon Detector providing Homology-based, Structural De Novo,
    and consensus detection with configurable in-between spacing for FASTA discovery.
    """

    def __init__(self, min_between: int = 10, min_tir_length: int = 8):
        self.homology_detector = HomologyDetector()
        self.structural_detector = StructuralDetector(
            min_between=min_between,
            min_tir_length=min_tir_length
        )

    def detect_all(
        self,
        genome_seq: str,
        library: Optional[Dict[str, TransposableElement]] = None,
        min_between: Optional[int] = None,
        min_tir_length: Optional[int] = None
    ) -> Dict[str, List[DetectionResult]]:
        """
        Run both Homology and Structural detectors and return consolidated results.
        """
        results: Dict[str, List[DetectionResult]] = {
            "homology": [],
            "structural": [],
            "consensus": []
        }

        if library:
            results["homology"] = self.homology_detector.detect(genome_seq, library)

        results["structural"] = self.structural_detector.detect(
            genome_seq,
            min_between=min_between,
            min_tir_length=min_tir_length
        )

        consensus = []
        for h in results["homology"]:
            consensus.append(h)

        for s in results["structural"]:
            overlap = False
            for h in results["homology"]:
                overlap_len = max(0, min(s.end, h.end) - max(s.start, h.start) + 1)
                if overlap_len > 0.6 * s.length:
                    overlap = True
                    break
            if not overlap:
                consensus.append(s)

        results["consensus"] = consensus
        return results

    def detect_new_candidates_from_sequence(
        self,
        seq: str,
        min_between: int = 10,
        min_tir_length: int = 8
    ) -> List[DetectionResult]:
        """
        Dedicated de novo discovery for newly provided FASTA sequences.
        Ensures the region in between the TIRs has at least `min_between` bases (default: 10 bp).
        """
        return self.structural_detector.detect(
            seq,
            min_between=min_between,
            min_tir_length=min_tir_length
        )
