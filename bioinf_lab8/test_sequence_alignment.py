"""Unit tests for sequence alignment module (sequence_alignment.py)."""

import pytest
from sequence_alignment import needleman_wunsch, format_matrix_text


def test_identical_sequences():
    seq = "ATGC"
    res = needleman_wunsch(seq, seq, match=2, mismatch=-1, gap=-1)
    assert res["aligned_seq1"] == "ATGC"
    assert res["aligned_seq2"] == "ATGC"
    assert res["alignment_marks"] == "||||"
    assert res["score"] == 8
    assert res["stats"]["matches"] == 4
    assert res["stats"]["mismatches"] == 0
    assert res["stats"]["gaps"] == 0


def test_complete_mismatch():
    seq1 = "AAAA"
    seq2 = "CCCC"
    res = needleman_wunsch(seq1, seq2, match=1, mismatch=-1, gap=-2)
    assert res["stats"]["matches"] == 0
    assert res["stats"]["mismatches"] == 4
    assert res["score"] == -4


def test_alignment_with_gaps():
    seq1 = "GATTACA"
    seq2 = "GCATGCA"
    res = needleman_wunsch(seq1, seq2, match=1, mismatch=-1, gap=-1)
    assert res["score"] == 2
    assert res["stats"]["matches"] == 5
    assert res["stats"]["gaps"] == 2
    assert res["stats"]["mismatches"] == 1


def test_traceback_path_bounds():
    seq1 = "ACTG"
    seq2 = "AG"
    res = needleman_wunsch(seq1, seq2, match=1, mismatch=-1, gap=-1)
    path = res["path"]
    assert path[0] == (0, 0)
    assert path[-1] == (len(seq1), len(seq2))


def test_format_matrix_text():
    seq1 = "AG"
    seq2 = "A"
    res = needleman_wunsch(seq1, seq2, match=1, mismatch=-1, gap=-1)
    formatted = format_matrix_text(res["score_matrix"], seq1, seq2, res["path"])
    assert "A" in formatted
    assert "-" in formatted
