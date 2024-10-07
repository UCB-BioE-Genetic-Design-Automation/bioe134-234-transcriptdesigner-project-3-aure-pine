from genedesign.seq_utils.reverse_complement import reverse_complement

def hairpin_counter(sequence, min_stem=3, min_loop=4, max_loop=9):
    """
    Counts the number of potential hairpin structures in a DNA sequence and returns a simple linear
    representation of the hairpins (stem1(loop)stem2_rc), or None if no hairpins are found.

    Parameters:
        sequence (str): The DNA sequence to analyze.
        min_stem (int): Minimum number of bases in the stem for stable hairpin.
        min_loop (int): Minimum number of bases in the loop.
        max_loop (int): Maximum number of bases in the loop.

    Returns:
        tuple: (int, str or None)
            - The count of potential hairpin structures.
            - A single string showing the detected hairpins in the format 'stem1(loop)stem2_rc', or None if no hairpins are found.
    """
    count = 0
    seq_len = len(sequence)
    hairpin_string = ""

    # Iterate through each base to consider it as a start of a stem
    for i in range(seq_len):
        # Only consider end points that would fit within the loop constraints
        for j in range(i + min_stem + min_loop, min(i + min_stem + max_loop + 1, seq_len)):
            stem1 = sequence[i:i+min_stem]
            stem2 = sequence[j:j+min_stem]

            # Check if the stems are complementary (reverse complement match)
            stem2_rc = reverse_complement(stem2)

            if stem1 == stem2_rc:
                count += 1

                # Extract the loop sequence
                loop = sequence[i+min_stem:j]

                # Create the linear representation (now correctly reversed for output)
                hairpin_representation = f"{stem1}({loop}){stem2}"

                # Append the linear hairpin representation to the string
                hairpin_string += f"Hairpin {count}: {hairpin_representation}\n"

    # Return count and the formatted hairpin string, or None if no hairpins found
    return count, hairpin_string if count > 0 else None

def non_stupid_hairpin_counter(sequence, min_stem=3, min_loop=4, max_loop=9):
    """
    Counts the number of potential hairpin structures in a DNA sequence. Hairpins are common secondary structures
    in nucleic acids where a sequence of nucleotides can fold back on itself to form a double-stranded stem with a single-stranded loop.

    The algorithm searches for regions within the sequence where a segment can base-pair with its reverse complement separated by a loop.
    This function scans for such occurrences by examining every possible substring as a potential stem and ensuring the intervening
    sequence, which would form the loop, meets the specified length requirements.

    Parameters:
        sequence (str): The DNA sequence to analyze.
        min_stem (int): Minimum number of bases in the stem for stable hairpin. The stem must be able to form at least this many
                        complementary base pairs to be considered.
        min_loop (int): Minimum number of bases in the loop. This prevents the formation of overly tight hairpins which may not be biologically relevant.
        max_loop (int): Maximum number of bases in the loop. This constrains the loop to a realistic size, preventing unlikely structures.

    Returns:
        int: The count of potential hairpin structures detected, indicating regions where secondary structure formation might inhibit biological processes like transcription or translation.

    This method does not account for the thermodynamic stability of the predicted hairpins, focusing solely on their potential for formation based on sequence complementarity and specified geometrical constraints.
    """
    count = 0
    seq_len = len(sequence)

    # Iterate through each base to consider it as a start of a stem
    for i in range(seq_len):
        # Only consider end points that would fit within the loop constraints
        for j in range(i + min_stem + min_loop, min(i + min_stem + max_loop + 1, seq_len)):
            stem1 = sequence[i:i+min_stem]
            stem2 = sequence[j:j+min_stem]

            # Check if the stems are complementary (reverse complement match)
            stem2_rc = ''.join(['ATCG'['TAGC'.index(n)] for n in stem2[::-1]])

            if stem1 == stem2_rc:
                count += 1

    return count

def main():
    # Example usage
    count, hairpins = hairpin_counter("AAAAAAAAAAAAAAAAAAAAAAAAAAA")
    print("Zero hairpins:", count)
    if hairpins:
        print(hairpins)
    else:
        print("No hairpins found.")

    count, hairpins = hairpin_counter("AAAAACCCAAAAAAAAAAGGGAAAAAA")
    print("CCC-N10-GGG:", count)
    if hairpins:
        print(hairpins)
    else:
        print("No hairpins found.")

    count, hairpins = hairpin_counter("AAAAACCCAAAAAAAAAGGGAAAAAA")
    print("CCC-N9-GGG:", count)
    if hairpins:
        print(hairpins)
    else:
        print("No hairpins found.")

    count, hairpins = hairpin_counter("AAAAACCCCAAAAAAAAGGGGAAAAAA")
    print("CCCC-N8-GGGG:", count)
    if hairpins:
        print(hairpins)
    else:
        print("No hairpins found.")

    count, hairpins = hairpin_counter("AAAAACACGAAAAAAAACGTGAAAAAA")
    print("CACG-N8-CGTG:", count)
    if hairpins:
        print(hairpins)
    else:
        print("No hairpins found.")

    count, hairpins = hairpin_counter("AAAACCCCCAAAAAAAAGGGGGAAA")
    print("CCCCC-N8-GGGGG:", count)
    if hairpins:
        print(hairpins)
    else:
        print("No hairpins found.")

    count, hairpins = hairpin_counter("AAAACCCCCAAAAAAAGGGGGAAA")
    print("CCCCC-N7-GGGGG:", count)
    if hairpins:
        print(hairpins)
    else:
        print("No hairpins found.")


if __name__ == "__main__":
    main()
