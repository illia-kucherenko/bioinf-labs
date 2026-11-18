"""
Unit Test Suite
Bioinformatics Lab 6 - Transposable Element Simulation, Detection & Intersection Analysis
"""

import unittest
from transposon_simulator import (
    reverse_complement,
    random_dna,
    TransposableElement,
    TransposonSimulator
)
from transposon_detector import (
    TransposonDetector,
    HomologyDetector,
    StructuralDetector,
    find_maximal_exact_matches,
    smith_waterman
)
from evaluation import IntersectionAnalyzer


class TestTransposonSimulator(unittest.TestCase):
    """Test synthetic DNA and transposon generation."""

    def test_reverse_complement(self):
        seq = "ATGCGATC"
        self.assertEqual(reverse_complement(seq), "GATCGCAT")
        self.assertEqual(reverse_complement("AAAA"), "TTTT")
        self.assertEqual(reverse_complement(""), "")

    def test_random_dna(self):
        dna = random_dna(300, seed=123)
        self.assertEqual(len(dna), 300)
        self.assertTrue(all(b in "ACGT" for b in dna))

    def test_transposable_element_creation(self):
        te = TransposableElement.create("TE_Test", total_length=50, tir_length=8, seed=42)
        self.assertEqual(te.length, 50)
        self.assertEqual(len(te.left_tir), 8)
        self.assertEqual(len(te.right_tir), 8)
        self.assertEqual(te.right_tir, reverse_complement(te.left_tir))

    def test_library_generation(self):
        sim = TransposonSimulator(seed=10)
        lib3 = sim.generate_library(3)
        self.assertEqual(len(lib3), 3)

        lib4 = sim.generate_library(4)
        self.assertEqual(len(lib4), 4)
        for name, te in lib4.items():
            self.assertGreaterEqual(te.length, 35)
            self.assertLessEqual(te.length, 65)

    def test_scenario_non_intersecting(self):
        sim = TransposonSimulator(seed=100)
        res = sim.simulate_scenario_non_intersecting(genome_length=350, te_count=4)
        self.assertEqual(res.final_length, 350)
        self.assertEqual(len(res.ground_truth_intervals), 4)

        # Check intervals are strictly sorted and non-overlapping
        prev_end = -1
        for gt in res.ground_truth_intervals:
            self.assertGreater(gt.start, prev_end)
            self.assertGreaterEqual(gt.end, gt.start)
            prev_end = gt.end
            # Sequence matches actual genome
            self.assertEqual(res.genome_sequence[gt.start : gt.end + 1], gt.sequence)

    def test_scenario_2_nested_intersection(self):
        sim = TransposonSimulator(seed=200)
        res = sim.simulate_scenario_2_intersecting_nested(genome_length=350, te_count=3)
        self.assertEqual(res.final_length, 350)
        # Should have: TE1_5', TE2_full, TE1_3', TE3_full = 4 intervals
        self.assertEqual(len(res.ground_truth_intervals), 4)

        te1_5p = res.ground_truth_intervals[0]
        te2_full = res.ground_truth_intervals[1]
        te1_3p = res.ground_truth_intervals[2]

        self.assertEqual(te1_5p.element_name, "TE_Mariner_1")
        self.assertEqual(te1_5p.fragment_type, "5_prime_fragment")
        self.assertEqual(te2_full.element_name, "TE_PiggyBac_2")
        self.assertEqual(te2_full.fragment_type, "full")
        self.assertEqual(te1_3p.element_name, "TE_Mariner_1")
        self.assertEqual(te1_3p.fragment_type, "3_prime_fragment")

        # Check contiguous insertion: TE2 is strictly between TE1_5' and TE1_3'
        self.assertEqual(te2_full.start, te1_5p.end + 1)
        self.assertEqual(te1_3p.start, te2_full.end + 1)

    def test_scenario_2_overlapping_intersection(self):
        sim = TransposonSimulator(seed=250)
        res = sim.simulate_scenario_2_intersecting_overlapping(genome_length=350, te_count=3)
        self.assertEqual(res.final_length, 350)
        # TE1 prefix, TE2 full, TE3 full
        self.assertEqual(len(res.ground_truth_intervals), 3)
        out = TransposonDetector().detect_all(res.genome_sequence, res.transposon_library)
        bench = IntersectionAnalyzer.evaluate(res, out["homology"])
        self.assertEqual(bench.recall, 100.0)

    def test_scenario_3_nested_intersection(self):
        sim = TransposonSimulator(seed=300)
        res = sim.simulate_scenario_3_intersecting_nested(genome_length=350, te_count=4)
        self.assertEqual(res.final_length, 350)
        # 6 intervals: TE1_5', TE2_5', TE3_full, TE2_3', TE1_3', TE4_full
        self.assertEqual(len(res.ground_truth_intervals), 6)

        te1_5p, te2_5p, te3, te2_3p, te1_3p, te4 = res.ground_truth_intervals

        # Check nested nesting
        self.assertEqual(te2_5p.start, te1_5p.end + 1)
        self.assertEqual(te3.start, te2_5p.end + 1)
        self.assertEqual(te2_3p.start, te3.end + 1)
        self.assertEqual(te1_3p.start, te2_3p.end + 1)


class TestTransposonDetector(unittest.TestCase):
    """Test detection algorithms and accuracy."""

    def setUp(self):
        self.sim = TransposonSimulator(seed=42)
        self.detector = TransposonDetector()

    def test_maximal_exact_matches(self):
        genome = "AAACCCGGGTTT"
        query = "CCCGGG"
        matches = find_maximal_exact_matches(genome, query, min_len=5)
        self.assertEqual(len(matches), 1)
        g_s, g_e, q_s, q_e, mlen = matches[0]
        self.assertEqual((g_s, g_e), (3, 8))
        self.assertEqual((q_s, q_e), (0, 5))
        self.assertEqual(mlen, 6)

    def test_smith_waterman_local_alignment(self):
        target = "ACGTACGTGGATCCAGTTAAAGCTTCGATCGAACGTACGT"
        query = "GGATCCAGTTAAAGCTTCGATC"
        alns = smith_waterman(target, query)
        self.assertGreater(len(alns), 0)
        s1_s, s1_e, s2_s, s2_e, score, ident, a1, a2 = alns[0]
        self.assertEqual((s1_s, s1_e), (8, 29))
        self.assertEqual(ident, 100.0)

    def test_detection_non_intersecting(self):
        res = self.sim.simulate_scenario_non_intersecting(350, 4)
        out = self.detector.detect_all(res.genome_sequence, res.transposon_library)
        homology_dets = out["homology"]

        self.assertEqual(len(homology_dets), 4)
        for d in homology_dets:
            self.assertEqual(d.status, "INTACT")
            self.assertEqual(d.identity, 100.0)

        bench = IntersectionAnalyzer.evaluate(res, homology_dets)
        self.assertEqual(bench.recall, 100.0)
        self.assertEqual(bench.precision, 100.0)
        self.assertEqual(bench.mean_iou, 100.0)

    def test_detection_2_nested_intersection(self):
        res = self.sim.simulate_scenario_2_intersecting_nested(350, 3)
        out = self.detector.detect_all(res.genome_sequence, res.transposon_library)
        homology_dets = out["homology"]

        # Check chained locus was detected
        chained = [d for d in homology_dets if d.status == "CHAINED_INTERRUPTED"]
        self.assertEqual(len(chained), 1)
        self.assertIn("TE_PiggyBac_2", chained[0].interrupted_by)

        # Check all ground truth segments are recovered
        bench = IntersectionAnalyzer.evaluate(res, homology_dets)
        self.assertEqual(bench.recall, 100.0)
        self.assertEqual(bench.mean_iou, 100.0)

    def test_detection_3_nested_intersection(self):
        res = self.sim.simulate_scenario_3_intersecting_nested(350, 4)
        out = self.detector.detect_all(res.genome_sequence, res.transposon_library)
        homology_dets = out["homology"]

        bench = IntersectionAnalyzer.evaluate(res, homology_dets)
        # All 6 segments should be detected with 100% recall
        self.assertEqual(bench.recall, 100.0)
        self.assertGreaterEqual(bench.mean_iou, 95.0)

    def test_structural_de_novo_detection(self):
        res = self.sim.simulate_scenario_2_intersecting_nested(350, 3)
        struct_dets = self.detector.detect_all(res.genome_sequence)["structural"]

        # Structural detector should identify TIRs and nested structures
        self.assertGreaterEqual(len(struct_dets), 2)
        nested_parents = [d for d in struct_dets if d.status == "STRUCTURAL_NESTED_PARENT"]
        self.assertGreaterEqual(len(nested_parents), 1)


if __name__ == "__main__":
    unittest.main()
