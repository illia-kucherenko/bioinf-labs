"""Unit tests for DNA pattern generation module (dna_pattern.py)."""

import pytest
import os
from dna_pattern import (
    clean_sequence,
    calculate_cg_content,
    calculate_information_content,
    sliding_window_pattern,
    calculate_center_of_weight,
    load_fasta
)


def test_clean_sequence():
    raw = "  cggactgatctatct\n(A27)\r\nCGTAGCATCTA "
    cleaned = clean_sequence("CGGACTGATCTATCTAAAACGTAGCATCTA")
    assert cleaned == "CGGACTGATCTATCTAAAACGTAGCATCTA"
    assert clean_sequence("actg\nACTG 123") == "ACTGACTG"


def test_calculate_cg_content_normal():
    assert calculate_cg_content("AAAAAAAAAA") == 0.0
    assert calculate_cg_content("TTTTTTTTTT") == 0.0
    assert calculate_cg_content("GGGGCCCCCC") == 100.0
    assert calculate_cg_content("ACGT") == 50.0
    assert calculate_cg_content("") == 0.0

    seq15 = "CGGACTGATCTATCT"
    expected = (7 / 15) * 100.0
    assert pytest.approx(calculate_cg_content(seq15), rel=1e-5) == expected


def test_calculate_cg_content_sequence_sensitive():
    sub = "CGCG"
    total_cg = 30.0
    val = calculate_cg_content(sub, mode="sequence_sensitive", total_cg_pct=total_cg)
    assert val == (30.0 / 4) * 4

    with pytest.raises(ValueError):
        calculate_cg_content(sub, mode="sequence_sensitive", total_cg_pct=None)


def test_calculate_information_content_homopolymer():
    poly_a = "A" * 30
    ic_val = calculate_information_content(poly_a, metric="ic")
    assert pytest.approx(ic_val, rel=1e-5) == 100.0

    sigma_val = calculate_information_content(poly_a, metric="sigma")
    assert pytest.approx(sigma_val, rel=1e-5) == 0.0


def test_calculate_information_content_alternating():
    alt = "AT" * 15
    ic_val = calculate_information_content(alt, metric="ic")
    sigma_val = calculate_information_content(alt, metric="sigma")
    assert pytest.approx(ic_val + sigma_val, rel=1e-5) == 100.0
    assert 40.0 < ic_val < 60.0


def test_sliding_window_test_sequence():
    prefix = "CGGACTGATCTATCT"
    middle = "A" * 27
    suffix = "CGTAGCATCTATCGATCTATCTAGCGATCTA"
    s = prefix + middle + suffix
    assert len(s) == 73

    L = 30
    pattern_44 = sliding_window_pattern(s, window_size=L, step=1, metric="ic", include_endpoint=True)
    assert len(pattern_44) == 44

    pattern_43 = sliding_window_pattern(s, window_size=L, step=1, metric="ic", include_endpoint=False)
    assert len(pattern_43) == 43

    for pt in pattern_44:
        assert 0.0 <= pt["cg"] <= 100.0
        assert 0.0 <= pt["ic"] <= 100.0
        assert len(pt["window"]) == L

    cw_cg, cw_ic = calculate_center_of_weight(pattern_44)
    assert 0.0 <= cw_cg <= 100.0
    assert 0.0 <= cw_ic <= 100.0


def test_fasta_loader():
    fasta_path = os.path.join("data", "gene_promoters_complete.fa")
    if os.path.exists(fasta_path):
        promoters = load_fasta(fasta_path)
        assert len(promoters) > 0
        first_gene = next(iter(promoters))
        entry = promoters[first_gene]
        assert "sequence" in entry
        assert len(entry["sequence"]) > 0
