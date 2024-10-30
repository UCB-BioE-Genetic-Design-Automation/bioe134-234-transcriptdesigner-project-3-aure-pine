import pytest
import numpy as np
from genedesign.seq_utils.sample_codon import SampleCodon

# Define the fixture to initialize and load the codon usage data once
@pytest.fixture(scope="module")
def sampler():
    sampler = SampleCodon()
    sampler.initiate()  # This will read from the actual 'genedesign/data/codon_usage.txt'
    return sampler

def test_initiate_loads_data_correctly(sampler):
    """Test that the initiate method correctly loads codon usage data."""
    
    # Check that codon_probabilities is populated
    assert sampler.codon_probabilities, "codon_probabilities should not be empty after initiation"
    
    # Check that known amino acids have data
    amino_acid = 'A'
    codons, probabilities = sampler.get_data(amino_acid)
    assert len(codons) > 0, f"No codons found for amino acid '{amino_acid}'"
    assert len(probabilities) > 0, f"No probabilities found for amino acid '{amino_acid}'"
    assert np.isclose(np.sum(probabilities), 1.0), "Probabilities should sum to 1"

def test_run_returns_codon_from_correct_set(sampler):
    """Test that the run method returns a codon from the correct set for a given amino acid."""
    amino_acid = 'A'
    codons = sampler.get_codons(amino_acid)
    codon = sampler.run(amino_acid)
    assert codon in codons, f"Codon '{codon}' not in expected codon list for amino acid '{amino_acid}'"

def test_run_raises_value_error_for_invalid_amino_acid(sampler):
    """Test that the run method raises a ValueError when an invalid amino acid is provided."""
    with pytest.raises(ValueError):
        sampler.run('Z')  # 'Z' is not a valid amino acid code

def test_run_sampling_probabilities(sampler):
    """Test that the run method samples codons according to the specified probabilities."""
    amino_acid = 'S'
    codons, expected_probs = sampler.get_data(amino_acid)
    num_samples = 100000
    samples = [sampler.run(amino_acid) for _ in range(num_samples)]
    counts = {codon: samples.count(codon) for codon in codons}
    observed_probs = [counts[codon] / num_samples for codon in codons]
    
    # Compare observed probabilities with expected probabilities
    for observed, expected in zip(observed_probs, expected_probs):
        assert abs(observed - expected) < 0.01, f"Observed probability {observed} differs from expected {expected}"

def test_run_with_single_codon_amino_acid(sampler):
    """Test that the run method correctly handles amino acids with a single codon."""
    amino_acid = 'M'
    codon = sampler.run(amino_acid)
    assert codon == 'ATG', "Codon for amino acid 'M' should be 'ATG'"

def test_run_with_stop_codon(sampler):
    """Test that the run method correctly samples stop codons."""
    amino_acid = '*'
    codons = sampler.get_codons(amino_acid)
    codon = sampler.run(amino_acid)
    assert codon in codons, f"Codon '{codon}' not in expected stop codon list"

def test_datastructure_method(sampler):
    """Test that the datastructure method returns the correct codon probabilities dictionary."""
    codon_probs = sampler.datastructure()
    assert isinstance(codon_probs, dict), "datastructure should return a dictionary"
    assert 'A' in codon_probs, "Amino acid 'A' should be in the codon probabilities dictionary"

def test_get_codons_method(sampler):
    """Test that the get_codons method returns the correct list of codons for an amino acid."""
    amino_acid = 'A'
    codons = sampler.get_codons(amino_acid)
    assert len(codons) > 0, f"No codons found for amino acid '{amino_acid}'"

def test_get_usages_method(sampler):
    """Test that the get_usages method returns the correct probabilities for an amino acid."""
    amino_acid = 'A'
    usages = sampler.get_usages(amino_acid)
    assert len(usages) > 0, f"No usages found for amino acid '{amino_acid}'"
    assert np.isclose(np.sum(usages), 1.0), "Probabilities should sum to 1"

def test_get_data_method(sampler):
    """Test that the get_data method returns the correct codons and usages for an amino acid."""
    amino_acid = 'A'
    codons, usages = sampler.get_data(amino_acid)
    assert len(codons) > 0, f"No codons found for amino acid '{amino_acid}'"
    assert len(usages) > 0, f"No usages found for amino acid '{amino_acid}'"
    assert np.isclose(np.sum(usages), 1.0), "Probabilities should sum to 1"

def test_total_probability_sums_to_one(sampler):
    """Test that the probabilities for each amino acid sum to 1."""
    for amino_acid in sampler.amino_acids:
        codons, usages = sampler.get_data(amino_acid)
        if len(usages) > 0:
            total_prob = np.sum(usages)
            assert np.isclose(total_prob, 1.0), f"Total probability for amino acid '{amino_acid}' does not sum to 1"