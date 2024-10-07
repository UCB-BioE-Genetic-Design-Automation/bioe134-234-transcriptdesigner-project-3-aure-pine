import pytest
from genedesign.seq_utils.check_seq import CheckSequence


@pytest.fixture
def check_sequence():
    """
    Fixture to initialize the CheckSequence object with real checkers.
    """
    cs = CheckSequence()
    cs.initiate()
    return cs


def test_all_checkers_pass(check_sequence):
    """
    Test case where all checkers pass.
    Assumes the sequence ["ATG", "GGT", "TAA"] passes all checks.
    """
    codons = ["ATG", "GGT", "GCG", "TAA"]  # Replace with codons that pass all checks
    result, results_list = check_sequence.run(codons)

    # Assert that all checkers passed
    print(f"Test all_checkers_pass results: {results_list}")
    assert result == True
    


def test_forbidden_sequence_fails(check_sequence):
    """
    Test case where the forbidden sequence checker fails.
    Assumes the sequence ["ATG", "GAA", "TTC", "TAA"] fails the forbidden sequence checker.
    """
    codons = ["ATG", "GAA", "TTC", "TAA"]  # EcoRI check
    result, results_list = check_sequence.run(codons)

    # Assert that the result is False and print the results_list for debugging
    print(f"Test forbidden_sequence_fails results: {results_list}")
    assert result == False


def test_promoter_sequence_fails(check_sequence):
    """
    Test case where the internal promoter checker fails.
    Assumes the sequence ["PROMOTER", "GGT", "TAA"] fails the promoter checker.
    """
    codons = ['ATG', 'TTG', 'ACA', 'GCT', 'AGC', 'TCA', 'GTC', 'CTA', 'GGT', 'ATA', "GGT", "TAA"]  # Codons that trigger the promoter checker
    result, results_list = check_sequence.run(codons)

    # Assert that the result is False and print the results_list for debugging
    print(f"Test forbidden_sequence_fails results: {results_list}")
    assert result == False
    

def test_hairpin_structure_fails(check_sequence):
    """
    Test case where the hairpin checker fails.
    Assumes the sequence ["ATG", "GGT", "HAIRPIN"] fails the hairpin checker.
    'CCC' ,'CCT' ,'TTC' ,'CCC' ,'CCA' ,'AAC' ,'CCC'
    """
    codons = ["ATG", "GGT", 'CCC' ,'CCT' ,'TTC' ,'CCC' ,'CCA' ,'AAC' ,'CCC']  # Replace with codons that trigger the hairpin checker
    result, results_list = check_sequence.run(codons)

    # Assert that the result is False and print the results_list for debugging
    print(f"Test forbidden_sequence_fails results: {results_list}")
    assert result == False
    


def test_gc_content_fails(check_sequence):
    """
    Test case where the GC content checker fails.
    Assumes the sequence ["ATG", "GC-RICH", "TAA"] fails the GC content checker.
    """
    codons = ["ATG", "GCG", 'CGT', 'GCG', 'CGC', "TAA"]  # Replace with codons that cause the GC content to fail
    result, results_list = check_sequence.run(codons)

    # Assert that the result is False and print the results_list for debugging
    assert result == False
    print(f"Test forbidden_sequence_fails results: {results_list}")


def test_multiple_checkers_fail(check_sequence):
    """
    Test case where multiple checkers fail.
    Assumes the sequence ["FORBIDDEN", "PROMOTER", "HAIRPIN"] fails multiple checkers.
    """
    codons = ["ATG", "GAA", "TTC", 'TTG', 'ACA', 'GCT', 'AGC', 'TCA', 'GTC', 'CTA', 'GGT', 
              'ATA', "GGT", 'CCC' ,'CCT' ,'TTC' ,'CCC' ,'CCA' ,'AAC' ,'CCC', "TAA"]
    result, results_list = check_sequence.run(codons)

    # Assert that the result is False and print the results_list for debugging
    assert result == False
    print(f"Test forbidden_sequence_fails results: {results_list}")