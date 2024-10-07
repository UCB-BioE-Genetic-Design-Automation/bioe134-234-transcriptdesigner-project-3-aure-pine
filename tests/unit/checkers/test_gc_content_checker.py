import pytest
from genedesign.checkers.gc_content_checker import gc_content_checker

def test_typical_case():
    # Case with 50% GC content
    dna_sequence = "ATGCGCATCG"
    assert gc_content_checker(dna_sequence) == True

def test_low_gc_content():
    # Case with 0% GC content
    dna_sequence = "ATATATAT"
    assert gc_content_checker(dna_sequence) == False

def test_high_gc_content():
    # Case with 100% GC content
    dna_sequence = "GGGGGGGG"
    assert gc_content_checker(dna_sequence) == False

def test_edge_case_empty_sequence():
    # Edge case with empty string
    with pytest.raises(ValueError, match="The DNA sequence cannot be empty."):
        gc_content_checker("")

def test_lower_case():
    dna_sequence = "ATGCGCATCG".lower()
    assert gc_content_checker(dna_sequence) == True

def test_gc_content_lower_boundary():
    # Case with exactly 40% GC content
    dna_sequence = "ATGCGTAGAT"
    assert gc_content_checker(dna_sequence) == True

def test_gc_content_upper_boundary():
    # Case with exactly 60% GC content
    dna_sequence = "ATGCGCGCTA"
    assert gc_content_checker(dna_sequence) == True

def test_gc_content_below_40_percent():
    # Case with slightly below 40% GC content
    dna_sequence = "ATGCGTATCTTA"
    assert gc_content_checker(dna_sequence) == False

def test_gc_content_above_60_percent():
    # Case with slightly above 60% GC content
    dna_sequence = "ATGCGCGCGCGA"
    assert gc_content_checker(dna_sequence) == False