from genedesign.seq_utils.sample_codon import SampleCodon
from genedesign.seq_utils.check_seq import CheckSequence

def montecarlo(window: str, last_n_codons: list[str], codon_sampler: SampleCodon, seq_checker: CheckSequence, n_codons_in_scope=3) -> list[str]:
    """
    Description:
    Generates a sequence of codons based on a window of amino acids and checks if the generated codons form a valid sequence using a sequence checker. 
    The function repeatedly generates codons until a valid sequence is found.

    Input:
    - window (str): A string of amino acids to be translated into codons.
    - last_n_codons (list[str]): A list of previously generated codons, used to maintain context when generating new codons.
    - codon_sampler (SampleCodon): An instance of a codon sampler that samples codons based on the given amino acids.
    - seq_checker (CheckSequence): An instance of a sequence checker that validates the concatenated codons (previous codons + generated codons) and returns a tuple (bool, info).
    - n_codons_in_scope (int, optional): The number of codons to return from the generated list. Defaults to 3.

    Output:
    - list[str]: A list of the first `n_codons_in_scope` valid codons that were generated.
    """
    # Ensure the window is not empty
    if not window:
        raise ValueError("Window cannot be empty")

    good_seq = False
    generated_codons = []
    attempts = 0
    max_attempts = 1000  # Limit attempts to prevent infinite loop

    # Continue sampling codons and checking the sequence until a valid one is generated
    while not good_seq and attempts < max_attempts:
        generated_codons = [codon_sampler.run(amino_acid) for amino_acid in window]
        good_seq, _ = seq_checker.run(last_n_codons + generated_codons)
        attempts += 1

    if not good_seq:
        raise RuntimeError("Unable to generate a valid sequence after multiple attempts")
    
    # Return only the first n codons in scope
    return generated_codons[:n_codons_in_scope]