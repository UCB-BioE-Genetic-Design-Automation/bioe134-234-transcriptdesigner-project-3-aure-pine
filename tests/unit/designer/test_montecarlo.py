import pytest
from unittest.mock import MagicMock, patch

# Assuming the MonteCarlo class is in a module named monte_carlo_module
from genedesign.montecarlo import MonteCarlo
# from genedesign.seq_utils.sample_codon import SampleCodon   
# from genedesign.rbs_chooser import rb

# Mocking external dependencies
@pytest.fixture
def mock_sample_codon():
    with patch('genedesign.seq_utils.sample_codon.SampleCodon') as MockSampleCodon:
        instance = MockSampleCodon.return_value
        instance.initiate.return_value = None
        # Mock the 'run' method to return a dummy codon (e.g., 'ATG')
        instance.run.return_value = 'ATG'
        yield instance

@pytest.fixture
def mock_rbs_chooser():
    with patch('genedesign.rbs_chooser.RBSChooser') as MockRBSChooser:
        instance = MockRBSChooser.return_value
        instance.initiate.return_value = None
        # Mock the 'optimized_run' method to return a dummy RBSOption
        instance.optimized_run.return_value = MagicMock(name='RBSOption')
        yield instance

@pytest.fixture
def mock_check_sequence():
    with patch('genedesign.seq_utils.check_seq.CheckSequence') as MockCheckSequence:
        instance = MockCheckSequence.return_value
        instance.initiate.return_value = None
        # Mock the 'run' method to return (True, 1.0)
        instance.run.return_value = (True, 1.0)
        yield instance

@pytest.fixture
def mock_codon_checker():
    with patch('genedesign.checkers.codon_checker.CodonChecker') as MockCodonChecker:
        instance = MockCodonChecker.return_value
        instance.initiate.return_value = None
        yield instance

@pytest.fixture
def monte_carlo(mock_sample_codon, mock_rbs_chooser, mock_check_sequence, mock_codon_checker):
    # Initialize the MonteCarlo instance
    monte_carlo_instance = MonteCarlo()
    monte_carlo_instance.initiate()
    return monte_carlo_instance

def test_run_with_valid_peptide(monte_carlo):
    peptide = 'MKTIIALSYIFCLVFADYKDDDDK'  # Example peptide sequence
    ignores = set()

    selected_rbs, codons = monte_carlo.run(peptide, ignores)

    # Assertions
    assert selected_rbs is not None, "Selected RBS should not be None"
    assert isinstance(codons, list), "Codons should be a list"
    assert len(codons) > 0, "Codons list should not be empty"
    # Check that codons are strings
    assert all(isinstance(codon, str) for codon in codons), "All codons should be strings"

def test_run_with_empty_peptide(monte_carlo):
    peptide = ''
    ignores = set()

    with pytest.raises(ValueError, match="Peptide needs to be a non-empty string."):
        monte_carlo.run(peptide, ignores)

def test_run_with_invalid_peptide_characters(monte_carlo):
    peptide = 'INVALID*PEPTIDE'
    ignores = set()

    # Depending on how your code handles invalid amino acids,
    # you may need to adjust this test.
    with pytest.raises(Exception):
        monte_carlo.run(peptide, ignores)

def test_montecarlo_with_special_first_codon(monte_carlo):
    # Mock the last_n_codons to be empty to trigger the special codon logic
    last_n_codons = []
    window = 'MKT'  # 'M' should be replaced with 'ATG'

    # Mock the sampler's 'run' method to return specific codons
    monte_carlo.sampler.run.side_effect = ['AAA', 'GGG']  # For 'K' and 'T'

    generated_codons = monte_carlo._MonteCarlo__montecarlo(window, last_n_codons)

    assert generated_codons[0] == 'ATG', "First codon should be 'ATG' for methionine"
    assert generated_codons[1:] == ['AAA', 'GGG'], "Generated codons do not match expected values"

def test_find_codons_success(monte_carlo):
    window = 'KTI'
    codons = ['ATG']
    selected_rbs = MagicMock(name='RBSOption')
    len_peptide = 10

    # Mock the checker to return success on the first attempt
    monte_carlo.checker.run.return_value = (True, 1.0)

    window_codons = monte_carlo._MonteCarlo__find_codons(window, codons, selected_rbs, len_peptide)

    assert len(window_codons) == monte_carlo.n_codons_in_scope, "Incorrect number of codons generated"
    # Since the sampler is mocked to return 'ATG', all codons will be 'ATG'
    assert all(codon == 'ATG' for codon in window_codons), "Codons do not match expected values"

def test_find_codons_failure(monte_carlo):
    window = 'KTI'
    codons = ['ATG']
    selected_rbs = MagicMock(name='RBSOption')
    len_peptide = 10

    # Mock the checker to always return False
    monte_carlo.checker.run.return_value = (False, 0.5)

    with pytest.raises(Exception, match="No valid window sequence found after"):
        monte_carlo._MonteCarlo__find_codons(window, codons, selected_rbs, len_peptide)