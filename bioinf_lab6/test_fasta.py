"""
Unit Tests for FASTA Input/Output and De Novo Transposon Discovery
Bioinformatics Lab 6 - Verification of FASTA reading and in-between spacing constraint (>= 10 bp).
"""

import os
import unittest
from fasta_io import parse_fasta_string, read_fasta_file, list_fasta_folder
from transposon_detector import TransposonDetector, StructuralDetector


class TestFastaDetection(unittest.TestCase):
    """Test FASTA reading and candidate transposon detection with >= 10 bp in-between constraint."""

    def test_parse_fasta_string(self):
        content = """>test_seq A test sequence description
ACGTACGT
ACGTACGT
>second_seq
GGCCAATT
"""
        records = parse_fasta_string(content)
        self.assertEqual(len(records), 2)
        self.assertEqual(records[0].id, "test_seq")
        self.assertEqual(records[0].description, "A test sequence description")
        self.assertEqual(records[0].sequence, "ACGTACGTACGTACGT")
        self.assertEqual(records[0].length, 16)
        self.assertEqual(records[1].id, "second_seq")
        self.assertEqual(records[1].sequence, "GGCCAATT")

    def test_read_fasta_files_in_folder(self):
        files = list_fasta_folder("fasta")
        self.assertGreaterEqual(len(files), 2)
        sample1 = next(f for f in files if "sample1" in f)
        records = read_fasta_file(sample1)
        self.assertEqual(len(records), 1)
        self.assertGreater(records[0].length, 300)

    def test_min_between_constraint_satisfied(self):
        """
        Verify that candidate transposons detected with min_between=10
        strictly have at least 10 bases in-between the left and right TIRs.
        """
        records = read_fasta_file("fasta/sample2_minimal_10bp_spacer.fasta")
        self.assertEqual(len(records), 1)
        seq = records[0].sequence

        detector = TransposonDetector(min_between=10, min_tir_length=8)
        candidates = detector.detect_new_candidates_from_sequence(seq, min_between=10, min_tir_length=8)

        self.assertGreater(len(candidates), 0)
        for c in candidates:
            # Check internal sequence between TIRs is >= 10 bp
            self.assertGreaterEqual(
                c.in_between_len,
                10,
                f"Candidate {c.element_name} has {c.in_between_len} bp in-between, expected >= 10 bp"
            )
            # Check TIR complementarity
            left_tir = c.left_tir
            right_tir = c.right_tir
            from transposon_simulator import reverse_complement
            self.assertEqual(
                right_tir,
                reverse_complement(left_tir),
                f"TIRs do not match reverse complement for {c.element_name}"
            )

    def test_exact_10bp_spacer_element_detected(self):
        """
        Specifically construct an artificial sequence containing a known transposon
        with exactly 10 bases between TIRs, and assert it is detected.
        """
        # Flank + TIR_L (8bp) + Spacer (10bp) + TIR_R (8bp) + Flank
        tir_l = "ACAGTGCT"
        from transposon_simulator import reverse_complement
        tir_r = reverse_complement(tir_l)  # "AGCACTGT"
        spacer_10bp = "NNNNACGTNN"  # 10 bases
        te_seq = tir_l + spacer_10bp + tir_r  # 26 bp total

        host_flank_left = "GGGGCCCCAAAATTTT"
        host_flank_right = "TTTTAAAACCCCGGGG"
        full_seq = host_flank_left + te_seq + host_flank_right

        sd = StructuralDetector(min_between=10, min_tir_length=8)
        dets = sd.detect(full_seq, min_between=10, min_tir_length=8)

        # Expected element starts at len(host_flank_left) = 16, ends at 16 + 26 - 1 = 41
        found = any(
            d.start == 16 and d.end == 41 and d.in_between_len == 10
            for d in dets
        )
        self.assertTrue(found, "Failed to detect transposon with exact 10-base spacer.")


if __name__ == "__main__":
    unittest.main()
