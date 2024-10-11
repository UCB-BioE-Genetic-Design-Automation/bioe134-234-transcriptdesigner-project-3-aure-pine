import pytest
import numpy as np
from unittest.mock import mock_open, patch
from genedesign.seq_utils.sample_codon import SampleCodon

np.random.seed(seed=42)

# TESTS ON MOCKS
# Test that `run` samples codons correctly based on the provided CAI probabilities
def test_run_sampling():
    # Mock the initiate method to avoid file operations
    sampler = SampleCodon()
    sampler.codon_probabilities = {
        'A': (np.array(['GCT', 'GCC', 'GCA', 'GCG']), np.array([0.4, 0.3, 0.2, 0.1])),
        'C': (np.array(['TGC', 'TGT']), np.array([0.6, 0.4]))
    }

    # Use a fixed seed to make the test deterministic
    np.random.seed(42)
    
    # Test sampling for amino acid 'A'
    codon = sampler.run('A')
    assert codon in ['GCT', 'GCC', 'GCA', 'GCG']

    # Test sampling for amino acid 'C'
    codon = sampler.run('C')
    assert codon in ['TGC', 'TGT']

# Test invalid amino acid input handling
def test_run_invalid_amino_acid():
    sampler = SampleCodon()
    sampler.codon_probabilities = {
        'A': (np.array(['GCT', 'GCC', 'GCA', 'GCG']), np.array([0.4, 0.3, 0.2, 0.1])),
        'C': (np.array(['TGC', 'TGT']), np.array([0.6, 0.4]))
    }

    with pytest.raises(ValueError, match="Invalid amino acid: X."):
        sampler.run('X')  # Invalid amino acid 'X'

# Test that probabilities sum up to 1
def test_probability_sum():
    sampler = SampleCodon()
    sampler.codon_probabilities = {
        'A': (np.array(['GCT', 'GCC', 'GCA', 'GCG']), np.array([0.4, 0.3, 0.2, 0.1])),
        'C': (np.array(['TGC', 'TGT']), np.array([0.6, 0.4]))
    }
    
    # Test that probabilities for 'A' sum to 1
    assert np.isclose(np.sum(sampler.codon_probabilities['A'][1]), 1.0)

    # Test that probabilities for 'C' sum to 1
    assert np.isclose(np.sum(sampler.codon_probabilities['C'][1]), 1.0)

# Test that the codon data structure is returned correctly
def test_datastructure():
    sampler = SampleCodon()
    sampler.codon_probabilities = {
        'A': (np.array(['GCT', 'GCC', 'GCA', 'GCG']), np.array([0.4, 0.3, 0.2, 0.1])),
        'C': (np.array(['TGC', 'TGT']), np.array([0.6, 0.4]))
    }

    data_structure = sampler.datastructure()
    assert 'A' in data_structure
    assert 'C' in data_structure
    assert np.array_equal(data_structure['A'][0], np.array(['GCT', 'GCC', 'GCA', 'GCG']))
    assert np.array_equal(data_structure['C'][0], np.array(['TGC', 'TGT']))


"""
Running on the actual data
"""
# Define the fixture to initialize and load the codon usage data once
@pytest.fixture(scope="module")
def sample_codon():
    sampler = SampleCodon()
    sampler.initiate()  # This will read from the actual 'genedesign/data/codon_usage.txt'
    return sampler

# Test the initiate method using the actual file
def test_initiate_with_actual_file(sample_codon):
    # Example checks: Ensure that key amino acids have their codon data loaded correctly
    # You can expand this by verifying all amino acids in your actual data
    
    # Check 'A' for Alanine
    codons, probabilities = sample_codon.codon_probabilities['A']
    assert len(codons) > 0  # Ensure codons for 'A' were loaded
    assert np.isclose(np.sum(probabilities), 1.0)  # Ensure probabilities sum to 1
    
    # Check 'C' for Cysteine
    codons, probabilities = sample_codon.codon_probabilities['C']
    assert len(codons) > 0  # Ensure codons for 'C' were loaded
    assert np.isclose(np.sum(probabilities), 1.0)  # Ensure probabilities sum to 1

# Test the run method with the actual file and data
def test_run_with_actual_file(sample_codon):
    # Use a fixed seed to make the test deterministic
    np.random.seed(42)

    # Test sampling for Alanine ('A')
    codon = sample_codon.run('A')
    assert codon in sample_codon.codon_probabilities['A'][0]  # Ensure valid codon for 'A'

    # Test sampling for Cysteine ('C')
    codon = sample_codon.run('C')
    assert codon in sample_codon.codon_probabilities['C'][0]  # Ensure valid codon for 'C'

# Test that codon probabilities sum up to 1 for all amino acids
def test_probability_sum_with_actual_file(sample_codon):
    for amino_acid, (codons, probabilities) in sample_codon.codon_probabilities.items():
        if len(probabilities) > 0:  # Only check amino acids with codons
            assert np.isclose(np.sum(probabilities), 1.0), f"Probabilities for {amino_acid} do not sum to 1"