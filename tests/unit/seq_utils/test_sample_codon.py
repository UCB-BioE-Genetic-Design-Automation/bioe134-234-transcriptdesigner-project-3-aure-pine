import pytest
import numpy as np
from genedesign.seq_utils.sample_codon import SampleCodon

np.random.seed(seed=42)

@pytest.fixture
def sample_codon():
    """
    Fixture to initialize the SampleCodon with codon usage data.
    """
    c = SampleCodon()
    c.initiate()  # Load the codon usage data
    return c

def test_data_structure(sample_codon):
    structure = sample_codon.datastructure()
    
    met = 'M'
    codons, probs = structure[met]

    assert len(codons) == 1
    assert len(probs) == 1

    out = sample_codon.run(met)
    assert out == 'ATG'

    codons2, probs2 = structure['A']
    
    assert len(codons2) == 4
    assert len(probs2) == 4

    assert any(item in codons2 for item in ['GCT', 'GCC', 'GCA', 'GCG'])
    assert any(item in probs2 for item in [0.16, 0.27, 0.22, 0.35])

# def test_sampling(sample_codon):
