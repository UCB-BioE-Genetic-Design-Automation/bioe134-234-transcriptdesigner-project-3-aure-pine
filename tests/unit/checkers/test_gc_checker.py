import pytest
from genedesign.checkers.gc_content_checker import gc_checker

# Test cases
def test_gc_content_exact_40_percent():
    sequence = "GGCCATATAT"  # 4 GC out of 10 nucleotides
    in_range, gc_content = gc_checker(sequence)
    assert in_range == True
    assert gc_content == 0.4

def test_gc_content_exact_60_percent():
    sequence = "GGCCGGATAT"  # 6 GC out of 10 nucleotides
    in_range, gc_content = gc_checker(sequence)
    assert in_range == True
    assert gc_content == 0.6

def test_gc_content_50_percent():
    sequence = "AGCTAGCTAGCTAGCT"
    in_range, gc_content = gc_checker(sequence)
    assert in_range == True
    assert gc_content == 0.5

def test_empty_sequence():
    sequence = ""
    with pytest.raises(ValueError, match="The sequence is empty."):
        gc_checker(sequence)

def test_zero_gc_content():
    sequence = "ATATATAT"
    in_range, gc_content = gc_checker(sequence)
    assert in_range == False
    assert gc_content == 0.0

def test_hundred_percent_gc_content():
    sequence = "GGGGCCCC"
    in_range, gc_content = gc_checker(sequence)
    assert in_range == False
    assert gc_content == 1.0

def test_lowercase_sequence():
    sequence = "gcatgcatgcat"
    in_range, gc_content = gc_checker(sequence)
    assert in_range == True
    assert gc_content == 0.5

def test_sequence_with_invalid_characters():
    sequence = "AGCTXAGC"
    with pytest.raises(ValueError, match="Invalid characters in sequence:"):
        gc_checker(sequence)

def test_gc_content_below_40_percent():
    sequence = "AAATTTAATC"
    in_range, gc_content = gc_checker(sequence)
    assert in_range == False
    assert gc_content == 0.1

def test_gc_content_above_60_percent():
    sequence = "GGGCCCGGGT"
    in_range, gc_content = gc_checker(sequence)
    assert in_range == False
    assert gc_content == 0.9

def test_gc_content_with_only_invalid_characters():
    sequence = "XYZ123"
    with pytest.raises(ValueError, match="Invalid characters in sequence:"):
        gc_checker(sequence)

def test_gc_content_with_mixed_valid_and_invalid_characters():
    sequence = "ATGCXYZATGC"
    with pytest.raises(ValueError, match="Invalid characters in sequence:"):
        gc_checker(sequence)

def test_gc_content_with_whitespace():
    sequence = "  A T G C  "
    with pytest.raises(ValueError, match="Invalid characters in sequence:"):
        gc_checker(sequence)

def test_gc_content_with_special_characters():
    sequence = "ATGC!@#$%^&*()"
    with pytest.raises(ValueError, match="Invalid characters in sequence:"):
        gc_checker(sequence)

def test_gc_content_just_below_40_percent():
    # 399 GC out of 1000 nucleotides = 39.9%
    sequence = "G" * 200 + "C" * 199 + "A" * 300 + "T" * 301
    in_range, gc_content = gc_checker(sequence)
    assert in_range == False
    assert gc_content == pytest.approx(0.399, 0.0001)

def test_gc_content_just_above_40_percent():
    # 401 GC out of 1000 nucleotides = 40.1%
    sequence = "G" * 201 + "C" * 200 + "A" * 299 + "T" * 300
    in_range, gc_content = gc_checker(sequence)
    assert in_range == True
    assert gc_content == pytest.approx(0.401, 0.0001)

def test_gc_content_just_below_60_percent():
    # 599 GC out of 1000 nucleotides = 59.9%
    sequence = "G" * 300 + "C" * 299 + "A" * 200 + "T" * 201
    in_range, gc_content = gc_checker(sequence)
    assert in_range == True
    assert gc_content == pytest.approx(0.599, 0.0001)

def test_gc_content_just_above_60_percent():
    # 601 GC out of 1000 nucleotides = 60.1%
    sequence = "G" * 301 + "C" * 300 + "A" * 199 + "T" * 200
    in_range, gc_content = gc_checker(sequence)
    assert in_range == False
    assert gc_content == pytest.approx(0.601, 0.0001)