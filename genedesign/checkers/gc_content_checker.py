def gc_content_checker(dna_seq: str) -> bool:
    """
    Description:
    Checks if the GC content of a given DNA sequence is between 40% and 60%.

    Input:
    dna_seq (str): A string representing the DNA sequence to be checked.

    Output:
    bool: Returns True if GC content is between 40% and 60%, otherwise False.

    Tests:
    """
    # Error handling for empty DNA sequence
    

    if not dna_seq:
        raise ValueError("The DNA sequence cannot be empty.")

    # Convert the sequence to uppercase to handle both lower and uppercase inputs
    dna_seq = dna_seq.upper()

    # Count the number of G's and C's
    g_count = dna_seq.count('G')
    c_count = dna_seq.count('C')

    # Calculate the length of the DNA
    total_length = len(dna_seq)

    # Calculate GC content percentage
    gc_content = (g_count + c_count) / total_length * 100

    # Return True if GC content is between 40% and 60%
    return 40 <= gc_content <= 60

