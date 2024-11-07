import re

def gc_checker(sequence):
    sequence = sequence.upper()

    # Check if the sequence is empty
    if not sequence:
        raise ValueError("The sequence is empty.")

    # Check for invalid characters
    valid_nucleotides = {'A', 'T', 'G', 'C'}
    invalid_chars = set(sequence) - valid_nucleotides
    if invalid_chars:
        raise ValueError(f"Invalid characters in sequence: {invalid_chars}")
    
    # Calculate the total number of nucleotides
    total_nucleotides = len(sequence)
    
    # Handle the case where the sequence might be empty
    if total_nucleotides == 0:
        gc_content = 0.0
    else:
        # Count the number of 'G' and 'C' nucleotides
        g_count = sequence.count('G')
        c_count = sequence.count('C')
        gc_count = g_count + c_count
        
        # Calculate GC content as a percentage
        gc_content = (gc_count / total_nucleotides)

    lower_bound = .4
    upper_bound = .6
    
    # Determine if GC content is between 40% and 60% inclusive
    in_range = lower_bound <= gc_content <= upper_bound
    
    return in_range, gc_content

if __name__ == '__main__':
    # Sample DNA sequences for testing
    sequences = [
        "AGCTATAG",          # GC content: 37.5%
        "GGCC",              # GC content: 100%
        "ATATATAT",          # GC content: 0%
        "",                  # Empty sequence
        "AGCTAGCTAGCTAGCT",  # GC content: 50%
        "gcatgcatgc",        # Lowercase letters, GC content: 50%
        "AGCTXAGC",          # Contains invalid character 'X'
        "CCCCGGGG",          # GC content: 100%
        "AAAATTTT",          # GC content: 0%
        "ACGTACGT",          # GC content: 50%
        "GGCCATATAT",        # GC content: 40%
        "GGCCGGATAT",        # GC content: 60%
    ]
    
    for seq in sequences:
        in_range, gc_content = gc_checker(seq)
        gc_percent = gc_content * 100  # Convert to percentage
        
        print(f"Sequence: '{seq}'")
        print(f"GC Content: {gc_percent:.2f}%")
        print(f"Is GC content between 40% and 60%? {in_range}")
        print("-" * 50)