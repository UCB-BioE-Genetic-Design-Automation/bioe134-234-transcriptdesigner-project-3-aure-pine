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

    # Calculate the GC content as an integer ratio (to avoid floating-point precision issues)
    gc_content_percentage = (g_count + c_count) * 100 / total_length

    # Return True if GC content is between 40% and 60%
    return 40 <= round(gc_content_percentage, 2) <= 60

