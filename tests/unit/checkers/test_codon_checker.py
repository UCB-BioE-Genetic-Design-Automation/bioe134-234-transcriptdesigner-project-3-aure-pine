from collections import Counter
import pytest
from genedesign.checkers.codon_checker import CodonChecker

@pytest.fixture
def codon_checker():
    """
    Fixture to initialize the CodonChecker with codon usage data.
    """
    checker = CodonChecker()
    checker.initiate()  # Load the codon usage data
    return checker

def test_empty_cds(codon_checker):
    """
    Test both run and my_run methods with an empty CDS.
    """
    cds = []
    len_peptide = 0  # Since CDS is empty

    expected_output = (False, 0.0, 0, 0.0)

    run_output = codon_checker.run(cds)
    my_run_output = codon_checker.my_run(cds, len_peptide)

    assert run_output == expected_output, "run method failed on empty CDS"
    assert my_run_output == expected_output, "my_run method failed on empty CDS"

def test_all_common_codons(codon_checker):
    """
    Test with a CDS containing only common codons (frequency >= 0.1).
    """
    # Example common codons based on typical usage frequencies
    cds = ['ATG', 'AAC', 'GAC', 'TGC', 'TAC', 'CAC', 'TTC', 'ATC', 'AAG', 'GAG', 'CAG']
    len_peptide = len(cds)  # Assuming one peptide per codon

    # Calculate expected values based on actual codon frequencies
    codon_counts = Counter(cds)
    codon_diversity = len(codon_counts) / 62  # Assuming 62 different codons
    print(codon_counts)
    print(codon_diversity)

    rare_codon_count = sum(codon_counts[codon] for codon in codon_checker.rare_codons if codon in cds)

    # CAI calculation as geometric mean of codon frequencies
    cai_numerators = [codon_checker.codon_frequencies.get(codon, 0.01) for codon in cds]
    cai_product = 1
    for freq in cai_numerators:
        cai_product *= freq
    cai_value = cai_product ** (1 / len(cai_numerators)) if cai_numerators else 0.0

    # Determine if codons are above board based on thresholds
    codons_above_board = (
        codon_diversity >= 0.5 and
        rare_codon_count <= 3 and
        cai_value >= 0.2
    )

    expected_output = (codons_above_board, codon_diversity, rare_codon_count, cai_value)

    run_output = codon_checker.run(cds)
    my_run_output = codon_checker.my_run(cds, len_peptide)
    print(run_output)
    print(my_run_output)

    assert run_output == expected_output, "run method failed on all common codons CDS"
    assert my_run_output == expected_output, "my_run method failed on all common codons CDS"

def test_with_rare_codons(codon_checker):
    """
    Test with a CDS containing some rare codons.
    """
    # Example CDS with rare codons (frequency < 0.1)
    cds = ['ATG', 'AGG', 'AGA', 'CTA', 'ATA', 'CGA', 'CCC']
    len_peptide = len(cds)  # Assuming one peptide per codon

    codon_counts = Counter(cds)
    codon_diversity = len(codon_counts) / 62

    # Calculate rare codon count
    rare_codon_count = sum(codon_counts[codon] for codon in codon_checker.rare_codons if codon in cds)

    # CAI calculation
    cai_numerators = [codon_checker.codon_frequencies.get(codon, 0.01) for codon in cds]
    cai_product = 1
    for freq in cai_numerators:
        cai_product *= freq
    cai_value = cai_product ** (1 / len(cai_numerators)) if cai_numerators else 0.0

    # Determine if codons are above board
    codons_above_board = (
        codon_diversity >= 0.5 and
        rare_codon_count <= 3 and
        cai_value >= 0.2
    )

    expected_output = (codons_above_board, codon_diversity, rare_codon_count, cai_value)

    run_output = codon_checker.run(cds)
    my_run_output = codon_checker.my_run(cds, len_peptide)

    assert run_output == expected_output, "run method failed on CDS with rare codons"
    assert my_run_output == expected_output, "my_run method failed on CDS with rare codons"

def test_with_unknown_codons(codon_checker):
    """
    Test with a CDS containing codons not present in codon_frequencies.
    """
    cds = ['ATG', 'XXX', 'YYY', 'GGC', 'AAA']
    len_peptide = len(cds)  # Assuming one peptide per codon

    codon_counts = Counter(cds)
    codon_diversity = len(codon_counts) / 62

    # Assume 'XXX' and 'YYY' are not in rare_codons
    rare_codon_count = sum(codon_counts[codon] for codon in codon_checker.rare_codons if codon in cds)

    # CAI calculation with unknown codons treated as frequency 0.01
    cai_numerators = [codon_checker.codon_frequencies.get(codon, 0.01) for codon in cds]
    cai_product = 1
    for freq in cai_numerators:
        cai_product *= freq
    cai_value = cai_product ** (1 / len(cai_numerators)) if cai_numerators else 0.0

    # Determine if codons are above board
    codons_above_board = (
        codon_diversity >= 0.5 and
        rare_codon_count <= 3 and
        cai_value >= 0.2
    )

    expected_output = (codons_above_board, codon_diversity, rare_codon_count, cai_value)

    run_output = codon_checker.run(cds)
    my_run_output = codon_checker.my_run(cds, len_peptide)

    assert run_output == expected_output, "run method failed on CDS with unknown codons"
    assert my_run_output == expected_output, "my_run method failed on CDS with unknown codons"

def test_thresholds_edge_case(codon_checker):
    """
    Test CDS that is exactly on the threshold values.
    """
    # Create a CDS with diversity just below and just above 0.5
    # For example, 31 unique codons (31/62 = 0.5)
    # To achieve exactly 0.5, use 31 unique codons each appearing once
    # Here, we'll use 31 unique codons repeated appropriately
    unique_codons = list(codon_checker.codon_frequencies.keys())[:31]  # Take first 31 codons
    cds = unique_codons  # Each codon appears once
    len_peptide = len(cds)  # Assuming one peptide per codon

    codon_counts = Counter(cds)
    codon_diversity = len(codon_counts) / 62  # Should be 31/62 = 0.5

    # Assume none of the first 31 codons are rare
    rare_codon_count = sum(codon_counts[codon] for codon in codon_checker.rare_codons if codon in cds)

    # CAI calculation
    cai_numerators = [codon_checker.codon_frequencies.get(codon, 0.01) for codon in cds]
    cai_product = 1
    for freq in cai_numerators:
        cai_product *= freq
    cai_value = cai_product ** (1 / len(cai_numerators)) if cai_numerators else 0.0

    # Determine if codons are above board
    codons_above_board = (
        codon_diversity >= 0.5 and
        rare_codon_count <= 3 and
        cai_value >= 0.2
    )

    # Adjust codon_diversity to be exactly 0.5
    assert codon_diversity == 0.5, "Codon diversity is not exactly 0.5 for the edge case"

    expected_output = (codons_above_board, codon_diversity, rare_codon_count, cai_value)

    run_output = codon_checker.run(cds)
    my_run_output = codon_checker.my_run(cds, len_peptide)

    assert run_output == expected_output, "run method failed on threshold edge case CDS"
    assert my_run_output == expected_output, "my_run method failed on threshold edge case CDS"

def test_varied_cds(codon_checker):
    """
    Test multiple CDS sequences with varied characteristics.
    """
    test_cases = [
        {
            'cds': ['ATG', 'TAT', 'TTT'],
            'len_peptide': 3,
            'description': 'Short CDS with all common codons',
        },
        {
            'cds': ['ATG', 'CTA', 'ATA', 'CGA'],
            'len_peptide': 4,
            'description': 'CDS with multiple rare codons',
        },
        {
            'cds': ['XXX', 'YYY', 'ZZZ'],
            'len_peptide': 3,
            'description': 'CDS with all unknown codons',
        },
        {
            'cds': ['ATG'] + ['TAT'] * 10 + ['ATA'] * 2 + ['CGA'] * 5,
            'len_peptide': 18,
            'description': 'Long CDS with some rare codons',
        },
    ]

    for case in test_cases:
        cds = case['cds']
        len_peptide = case['len_peptide']
        description = case['description']

        codon_counts = Counter(cds)
        codon_diversity = len(codon_counts) / 62

        # Calculate rare_codon_count
        rare_codon_count = sum(codon_counts[codon] for codon in codon_checker.rare_codons if codon in cds)

        # CAI calculation
        cai_numerators = [codon_checker.codon_frequencies.get(codon, 0.01) for codon in cds]
        cai_product = 1
        for freq in cai_numerators:
            cai_product *= freq
        cai_value = cai_product ** (1 / len(cai_numerators)) if cai_numerators else 0.0

        # Determine if codons are above board
        codons_above_board = (
            codon_diversity >= 0.5 and
            rare_codon_count <= 3 and
            cai_value >= 0.2
        )

        expected_output = (codons_above_board, codon_diversity, rare_codon_count, cai_value)

        run_output = codon_checker.run(cds)
        my_run_output = codon_checker.my_run(cds, len_peptide)

        assert run_output == expected_output, f"run method failed on test case: {description}"
        assert my_run_output == expected_output, f"my_run method failed on test case: {description}"