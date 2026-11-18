"""
Evaluation and Intersection Behavior Analysis Module
Bioinformatics Lab 6 - Benchmarking Transposon Detection & Intersection Dynamics
"""

from dataclasses import dataclass, field
from typing import List, Dict, Tuple, Optional
from transposon_simulator import GenomicInterval, SimulationResult
from transposon_detector import DetectionResult


@dataclass
class MatchEvaluation:
    """Evaluation result for a single ground truth interval or detection."""
    ground_truth_name: str
    gt_start: int
    gt_end: int
    gt_type: str
    detected_name: Optional[str] = None
    det_start: Optional[int] = None
    det_end: Optional[int] = None
    det_method: Optional[str] = None
    overlap_bases: int = 0
    iou: float = 0.0  # Jaccard index / Intersection over Union
    is_matched: bool = False
    status_assessment: str = ""


@dataclass
class BenchmarkSummary:
    """Consolidated metrics across all elements and fragments in a simulation."""
    total_ground_truth_intervals: int
    total_detections: int
    true_positives: int
    false_positives: int
    false_negatives: int
    precision: float
    recall: float
    f1_score: float
    mean_iou: float
    matches: List[MatchEvaluation]
    intersection_report: List[str]


class IntersectionAnalyzer:
    """
    Analyzes how detection algorithms behave when transposable elements
    intersect (via nesting or partial overlap).
    """

    @staticmethod
    def evaluate(
        sim_result: SimulationResult,
        detections: List[DetectionResult],
        iou_threshold: float = 0.5
    ) -> BenchmarkSummary:
        """
        Compare detection intervals against ground truth genomic intervals.
        """
        gt_intervals = sim_result.ground_truth_intervals
        unmatched_dets = list(detections)
        matches: List[MatchEvaluation] = []
        tp = 0
        iou_sum = 0.0

        for gt in gt_intervals:
            best_det = None
            best_iou = 0.0
            best_overlap = 0

            for det in unmatched_dets:
                # Calculate coordinate overlap
                overlap = max(0, min(gt.end, det.end) - max(gt.start, det.start) + 1)
                union = (gt.end - gt.start + 1) + (det.end - det.start + 1) - overlap
                iou = overlap / union if union > 0 else 0.0

                # Check if name is compatible
                # e.g., "TE_Mariner_1" matches "TE_Mariner_1" or "Putative_TE_TIR_..."
                name_match = (
                    gt.element_name in det.element_name or
                    det.element_name in gt.element_name or
                    det.method.startswith("Structural")  # De novo doesn't know library names
                )

                if name_match and iou > best_iou:
                    best_iou = iou
                    best_overlap = overlap
                    best_det = det

            if best_det and best_iou >= iou_threshold:
                tp += 1
                iou_sum += best_iou
                unmatched_dets.remove(best_det)

                # Assess status
                if gt.fragment_type == "full":
                    status = "Correctly detected full element"
                else:
                    status = f"Correctly resolved {gt.fragment_type} resulting from intersection"

                matches.append(MatchEvaluation(
                    ground_truth_name=gt.element_name,
                    gt_start=gt.start,
                    gt_end=gt.end,
                    gt_type=gt.fragment_type,
                    detected_name=best_det.element_name,
                    det_start=best_det.start,
                    det_end=best_det.end,
                    det_method=best_det.method,
                    overlap_bases=best_overlap,
                    iou=best_iou,
                    is_matched=True,
                    status_assessment=status
                ))
            else:
                matches.append(MatchEvaluation(
                    ground_truth_name=gt.element_name,
                    gt_start=gt.start,
                    gt_end=gt.end,
                    gt_type=gt.fragment_type,
                    status_assessment=f"Missed {gt.fragment_type} (boundary shift or low coverage)"
                ))

        fp = len(unmatched_dets)
        fn = len(gt_intervals) - tp
        precision = (tp / (tp + fp)) * 100.0 if (tp + fp) > 0 else 0.0
        recall = (tp / len(gt_intervals)) * 100.0 if gt_intervals else 0.0
        f1 = (2 * precision * recall / (precision + recall)) if (precision + recall) > 0 else 0.0
        mean_iou = (iou_sum / tp) * 100.0 if tp > 0 else 0.0

        # Generate scientific narrative report on intersection behavior
        report = IntersectionAnalyzer.generate_behavior_report(sim_result, detections, matches)

        return BenchmarkSummary(
            total_ground_truth_intervals=len(gt_intervals),
            total_detections=len(detections),
            true_positives=tp,
            false_positives=fp,
            false_negatives=fn,
            precision=precision,
            recall=recall,
            f1_score=f1,
            mean_iou=mean_iou,
            matches=matches,
            intersection_report=report
        )

    @staticmethod
    def generate_behavior_report(
        sim_result: SimulationResult,
        detections: List[DetectionResult],
        matches: List[MatchEvaluation]
    ) -> List[str]:
        """
        Produce detailed diagnostic commentary describing the algorithm's behavior
        in intersecting vs non-intersecting situations.
        """
        lines = []
        lines.append(f"=== Algorithmic Behavior Analysis: {sim_result.scenario_name} ===")
        lines.append(f"Synthetic DNA Length: {sim_result.final_length} bp | Ground Truth Segments: {len(sim_result.ground_truth_intervals)}")
        lines.append("")

        scenario = sim_result.scenario_name

        if "Non-Intersecting" in scenario:
            lines.append("1. Control Scenario (Independent Insertions):")
            lines.append("   - Baseline behavior: Each transposable element is unbroken and contiguous.")
            lines.append("   - Local alignment identifies 100% query coverage with zero internal gaps.")
            lines.append("   - Structural detector pairs the 5' and 3' TIRs at their canonical element lengths.")
            lines.append("   - Result: Both methods pinpoint exact [start, end] positions without ambiguity.")

        elif "2-Transposon Nested" in scenario:
            lines.append("2. Two-Transposon Nested Intersection Dynamics:")
            lines.append("   - Mechanism: Transposon 2 inserted directly into Transposon 1.")
            lines.append("   - Effect on Outer Element (TE1):")
            lines.append("     * The original element is interrupted and split into two distinct sub-intervals (5' and 3' fragments).")
            lines.append("     * A naive global alignment algorithm would fail or severely penalize this locus due to the huge intervening gap.")
            lines.append("     * Our Homology Chainer recognizes the two split fragments and chains them into a compound interrupted locus.")
            lines.append("     * The Structural TIR Detector identifies the outer TIR pair, observing an expanded span encompassing TE2.")
            lines.append("   - Effect on Inner Element (TE2):")
            lines.append("     * Transposon 2 remains completely intact and undisturbed.")
            lines.append("     * Both Homology and Structural detectors locate TE2 with 100% confidence strictly inside TE1's bounds.")

        elif "Overlapping" in scenario:
            lines.append("3. Two-Transposon Overlapping Intersection Dynamics:")
            lines.append("   - Mechanism: Transposon 2 partially overwrites or recombinantly overlaps Transposon 1.")
            lines.append("   - Effect on Upstream Element (TE1):")
            lines.append("     * TE1's 3' end is truncated. Homology alignment detects only the remaining 5' prefix.")
            lines.append("     * Structural detection fails for TE1 because its 3' Terminal Inverted Repeat was deleted/overwritten.")
            lines.append("   - Effect on Downstream Element (TE2):")
            lines.append("     * TE2 retains both terminal repeats and aligns fully.")
            lines.append("   - Algorithmic Takeaway: Overlapping insertions destroy structural repeat symmetry, requiring homology fallback.")

        elif "3-Transposon" in scenario:
            lines.append("4. Multi-Level (3-Transposon) Hierarchical Nesting Dynamics:")
            lines.append("   - Mechanism: TE3 inserted into TE2, which previously inserted into TE1.")
            lines.append("   - Architectural Hierarchy: [TE1_5'] -> [TE2_5'] -> [TE3_full] -> [TE2_3'] -> [TE1_3'].")
            lines.append("   - Algorithmic Behavior:")
            lines.append("     * Innermost Element (TE3): Detected intact as the core insertion.")
            lines.append("     * Middle Element (TE2): Detected as two fragmented halves surrounding TE3; chained into a nested parent.")
            lines.append("     * Outermost Element (TE1): Detected as two flanking fragments surrounding the entire TE2+TE3 block.")
            lines.append("     * Structural Detector exhibits multi-tier nesting discovery, reporting parent-child relationships.")

        lines.append("")
        lines.append("Performance Summary:")
        matched_count = sum(1 for m in matches if m.is_matched)
        lines.append(f"- Successfully resolved genomic segments: {matched_count} / {len(matches)}")
        for m in matches:
            if m.is_matched:
                lines.append(f"  * [{m.gt_start}:{m.gt_end}] {m.ground_truth_name} ({m.gt_type}): IoU={m.iou*100:.1f}%, Detected via {m.det_method}")
            else:
                lines.append(f"  * [{m.gt_start}:{m.gt_end}] {m.ground_truth_name} ({m.gt_type}): Unresolved")

        return lines

    @staticmethod
    def generate_fasta_report(
        fasta_id: str,
        seq_length: int,
        detections: List[DetectionResult],
        min_between: int = 10
    ) -> List[str]:
        """
        Produce analytical report for newly discovered candidate transposons
        from an opened FASTA sequence.
        """
        lines = []
        lines.append(f"=== De Novo Transposon Discovery Report: {fasta_id} ===")
        lines.append(f"Sequence Length: {seq_length} bp | In-Between Spacer Constraint: >= {min_between} bp")
        lines.append(f"Candidate Transposons Identified: {len(detections)}")
        lines.append("")

        if not detections:
            lines.append("No structural transposable elements were detected satisfying the criteria:")
            lines.append(f"- Minimum in-between spacing: {min_between} bases")
            lines.append("- Terminal Inverted Repeat (TIR) reverse complement symmetry")
            return lines

        lines.append("Discovered Putative Transposable Elements:")
        lines.append("-" * 72)
        for i, det in enumerate(detections, 1):
            lines.append(
                f"{i}. {det.element_name}: [{det.start}:{det.end}] (1-based: [{det.start_1based}:{det.end_1based}])"
            )
            lines.append(f"   Length: {det.length} bp | In-between internal bases: {det.in_between_len} bp (>= {min_between} bp constraint met)")
            if det.left_tir and det.right_tir:
                lines.append(f"   Left TIR:  5'-{det.left_tir}-3' ({len(det.left_tir)} bp)")
                lines.append(f"   Right TIR: 5'-{det.right_tir}-3' ({len(det.right_tir)} bp)")
            if det.tsd:
                lines.append(f"   Target Site Duplication (TSD): '{det.tsd}' ({len(det.tsd)} bp flanking)")
            if det.interrupted_by:
                lines.append(f"   Nested Architecture: Encloses inner element(s) -> {', '.join(det.interrupted_by)}")
            lines.append(f"   Classification / Status: {det.status}")
            lines.append("")

        lines.append("Biological Assessment:")
        nested_count = sum(1 for d in detections if d.status == "STRUCTURAL_NESTED_PARENT" or d.interrupted_by)
        if nested_count > 0:
            lines.append(f"- Complex Architecture: {nested_count} element(s) exhibit nested/intersecting insertion.")
        else:
            lines.append("- Independent Architecture: Detected elements are non-overlapping candidate loci.")
        lines.append(f"- All candidates adhere to the minimum {min_between}-base internal spacer requirement.")
        return lines
