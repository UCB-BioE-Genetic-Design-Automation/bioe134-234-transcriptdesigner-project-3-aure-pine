import pytest
from genedesign.sliding_window_generator import sliding_window_generator  # Replace 'your_module' with the actual module name

def test_sliding_window_typical_case():
    # Test a typical case with default parameters
    sequence = "MTKLPQACDGHI"
    result = list(sliding_window_generator(sequence))  # Convert generator to list

    # Check the expected output
    expected = [
        "MTKLPQACD",  # First window: 3 in-scope ("MTK"), 6 context-after ("LPQACD")
        "LPQACDGHI",  # Second window: 3 in-scope ("LPQ"), 6 context-after ("ACDGHI")
        "ACDGHI",     # Third window: 3 in-scope ("ACD"), only 3 context-after left ("GHI")
        "GHI"         # Last window: 3 in-scope ("GHI"), none left after.
    ]
    assert result == expected

def test_sliding_window_custom_step():
    # Test case with a custom step size
    sequence = "MTKLPQACDGHI"
    result = list(sliding_window_generator(sequence, step=1))  # Step size of 1

    expected = [
        "MTKLPQACD",  # First window
        "TKLPQACDG",  # Second window
        "KLPQACDGH",  # Third window
        "LPQACDGHI",  # Fourth window
        "PQACDGHI",   # Fifth window
        "QACDGHI",    # Sixth window
        "ACDGHI",     # Seventh window
        "CDGHI",      # Eighth window
        "DGHI",       # Ninth window
        "GHI",        # Tenth window
        "HI",         # Eleventh window
        "I"           # Last window (remaining single amino acid)
    ]
    assert result == expected

def test_sliding_window_small_sequence():
    # Test case where the sequence is smaller than the window
    sequence = "MTK"
    result = list(sliding_window_generator(sequence))  # Convert generator to list

    # Check the expected output (only one window, no context-after)
    expected = ["MTK"]
    assert result == expected

def test_sliding_window_large_n_in_scope():
    # Test case where n_in_scope is larger than the sequence
    sequence = "MTK"
    result = list(sliding_window_generator(sequence, n_in_scope=5))  # n_in_scope larger than sequence

    # Check the expected output (the entire sequence is the only window)
    expected = ["MTK"]
    assert result == expected

def test_sliding_window_large_n_ahead():
    # Test case with n_ahead = 10, and the sequence length is at least double n_ahead
    sequence = "MTKLPQACDGHIJKLMNOQRST"  # Sequence with 22 characters
    result = list(sliding_window_generator(sequence, n_ahead=10))  # n_ahead is 10

    # Check the expected output
    expected = [
        "MTKLPQACDGHIJ",   # First window: 3 in-scope ("MTK"), 10 context-after ("LPQACDGHIJ")
        "LPQACDGHIJKLM",  # Second window: 3 in-scope ("LPQ"), 10 context-after ("ACDGHIJKLM")
        "ACDGHIJKLMNOQ",     # Third window: 3 in-scope ("ACD"), 10 context-after ("GHIJKLMNOQ")
        "GHIJKLMNOQRST",        # Fourth window: 3 in-scope ("GHI"), remaining 10 context-after ("JKLMNOQRST")
        "JKLMNOQRST",           # Fifth window: 3 in-scope ("JKL"), remaining 7 context-after ("MNOQRST")
        "MNOQRST",              # Sixth window: 3 in-scope ("MNO"), remaining 4 context-after ("QRST")
        "QRST",                 # Seventh window: 3 in-scope ("QRS"), remaining 1 context-after ("T")
        "T"                     # Eighth window: 3 in-scope ("T"), no context-after remaining
    ]
    
    assert result == expected

def test_sliding_window_empty_sequence():
    # Test case with an empty sequence
    sequence = ""
    result = list(sliding_window_generator(sequence))  # Convert generator to list

    # Check the expected output (no windows)
    assert result == []

def test_sliding_window_custom_n_in_scope_n_ahead():
    # Test case with custom n_in_scope and n_ahead
    sequence = "MTKLPQACDGHI"
    result = list(sliding_window_generator(sequence, n_in_scope=4, n_ahead=2))  # Custom window sizes

    expected = [
        "MTKLPQ",  # First window: 4 in-scope ("MTKL"), 2 context-after ("PQ")
        "PQACDG",  # Second window: 4 in-scope ("PQAC"), 2 context-after ("DG")
        "DGHI"     # Third window: 4 in-scope ("DGHI"), no context-after left
    ]
    assert result == expected