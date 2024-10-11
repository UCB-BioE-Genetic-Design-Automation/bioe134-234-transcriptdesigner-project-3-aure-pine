import pytest
from unittest.mock import Mock
from genedesign.seq_utils.sample_codon import SampleCodon
from genedesign.seq_utils.check_seq import CheckSequence
from genedesign.montecarlo import montecarlo  # Replace 'your_module' with the actual module name

def test_montecarlo_typical_case():
    # Mock the SampleCodon instance
    codon_sampler = Mock(SampleCodon)
    codon_sampler.run.side_effect = ['ATG', 'ACA', 'AAA']  # Codons for "MTK"

    # Mock the CheckSequence instance
    seq_checker = Mock(CheckSequence)
    seq_checker.return_value = (True, 'Valid')  # Always return valid sequence

    # Test inputs
    window = "MTK"
    last_n_codons = ["ATG", "ACA"]
    n_codons_in_scope = 3

    # Run the function
    result = montecarlo(window, last_n_codons, codon_sampler, seq_checker, n_codons_in_scope)

    # Check the result
    assert result == ["ATG", "ACA", "AAA"]

def test_montecarlo_single_codon():
    # Mock the SampleCodon instance
    codon_sampler = Mock(SampleCodon)
    codon_sampler.run.side_effect = ['ATG']  # Codon for "M"

    # Mock the CheckSequence instance
    seq_checker = Mock(CheckSequence)
    seq_checker.return_value = (True, 'Valid')

    # Test inputs
    window = "M"
    last_n_codons = ["ATG"]
    n_codons_in_scope = 1

    # Run the function
    result = montecarlo(window, last_n_codons, codon_sampler, seq_checker, n_codons_in_scope)

    # Check the result
    assert result == ["ATG"]

def test_montecarlo_empty_window():
    # Mock the SampleCodon instance
    codon_sampler = Mock(SampleCodon)
    
    # Mock the CheckSequence instance
    seq_checker = Mock(CheckSequence)

    # Test inputs
    window = ""
    last_n_codons = ["ATG"]

    # Expect a ValueError for an empty window
    with pytest.raises(ValueError, match="Window cannot be empty"):
        montecarlo(window, last_n_codons, codon_sampler, seq_checker)

def test_montecarlo_max_attempts():
    # Mock the SampleCodon instance
    codon_sampler = Mock(SampleCodon)
    # Use a lambda to continually return the same codons, avoiding StopIteration
    codon_sampler.run.side_effect = lambda amino_acid: 'ATG'

    # Mock the CheckSequence instance to always return False, so it keeps looping
    seq_checker = Mock(CheckSequence)
    seq_checker.return_value = (False, 'Invalid')

    # Test inputs
    window = "MTK"
    last_n_codons = ["ATG", "ACA"]

    # Expect a RuntimeError because the sequence never becomes valid
    with pytest.raises(RuntimeError, match="Unable to generate a valid sequence after multiple attempts"):
        montecarlo(window, last_n_codons, codon_sampler, seq_checker)

def test_montecarlo_partial_codons_returned():
    # Mock the SampleCodon instance
    codon_sampler = Mock(SampleCodon)
    codon_sampler.run.side_effect = ['ATG', 'ACA', 'AAA']

    # Mock the CheckSequence instance to return True after the first two codons
    seq_checker = Mock(CheckSequence)
    seq_checker.return_value = (True, 'Valid')

    # Test inputs
    window = "MTK"
    last_n_codons = ["ATG", "ACA"]
    n_codons_in_scope = 2  # We only want the first 2 codons

    # Run the function
    result = montecarlo(window, last_n_codons, codon_sampler, seq_checker, n_codons_in_scope)

    # Check the result
    assert result == ["ATG", "ACA"]