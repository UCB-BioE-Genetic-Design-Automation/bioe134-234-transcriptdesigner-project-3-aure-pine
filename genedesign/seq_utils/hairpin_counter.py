from genedesign.seq_utils.reverse_complement import reverse_complement
import numpy as np
from numba import njit

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

def optimized_hairpin_counter(sequence, min_stem=3, min_loop=4, max_loop=9):
    """
    Counts the number of potential hairpin structures in a DNA sequence and returns a linear
    representation of all detected hairpins (stem1(loop)stem2_rc), or None if no hairpins are found.

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
    hairpin_string = []

    # Precompute reverse complements for all possible stems in the sequence
    reverse_complements = {}
    for i in range(seq_len - min_stem + 1):
        stem = sequence[i:i + min_stem]
        reverse_complements[stem] = reverse_complement(stem)

    # Iterate through each base to consider it as a start of a stem
    for i in range(seq_len - min_stem):
        stem1 = sequence[i:i + min_stem]

        # Only consider end points that would fit within the loop constraints
        for j in range(i + min_stem + min_loop, min(i + min_stem + max_loop + 1, seq_len - min_stem + 1)):
            stem2 = sequence[j:j + min_stem]

            # Use precomputed reverse complement of stem2
            stem2_rc = reverse_complements[stem2]

            # Check if the stems are complementary (reverse complement match)
            if stem1 == stem2_rc:
                count += 1

                # Extract the loop sequence
                loop = sequence[i + min_stem:j]

                # Create the linear representation (now correctly reversed for output)
                hairpin_representation = f"{stem1}({loop}){stem2}"

                # Append the linear hairpin representation to the list
                hairpin_string.append(f"Hairpin {count}: {hairpin_representation}")

                # Do not break here, continue to find more potential hairpins

    # Return count and the formatted hairpin string, or None if no hairpins found
    return count, '\n'.join(hairpin_string) if count > 0 else None

def optimized_non_stupid_hairpin_counter(sequence, min_stem=3, min_loop=4, max_loop=9):
    """
    Counts the number of potential hairpin structures in a DNA sequence using an optimized algorithm. Hairpins are common secondary structures
    in nucleic acids where a sequence of nucleotides can fold back on itself to form a double-stranded stem with a single-stranded loop.

    The optimized algorithm improves performance by representing nucleotides and stems as integers, precomputing reverse complements, and reducing redundant computations.
    It searches for regions within the sequence where a segment can base-pair with its reverse complement separated by a loop, scanning for such occurrences by examining
    every possible substring as a potential stem while ensuring the intervening sequence (the loop) meets the specified length requirements.

    Parameters:
        sequence (str): The DNA sequence to analyze.
        min_stem (int): Minimum number of bases in the stem for a stable hairpin. The stem must be able to form at least this many
                        complementary base pairs to be considered.
        min_loop (int): Minimum number of bases in the loop. This prevents the formation of overly tight hairpins which may not be biologically relevant.
        max_loop (int): Maximum number of bases in the loop. This constrains the loop to a realistic size, preventing unlikely structures.

    Returns:
        int: The count of potential hairpin structures detected, indicating regions where secondary structure formation might inhibit biological processes like transcription or translation.

    This method does not account for the thermodynamic stability of the predicted hairpins, focusing solely on their potential for formation based on sequence complementarity and specified geometrical constraints.

    Optimizations implemented:
        - Nucleotides are converted to integer representations for faster computation.
        - Reverse complements of all possible stems are precomputed to avoid redundant calculations.
        - Stem sequences are represented as integers to leverage efficient arithmetic operations and reduce memory usage.
    """
    count = 0
    seq_len = len(sequence)
    k = min_stem

    # Mapping nucleotides to integers
    base_to_int = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
    comp_int = [3, 2, 1, 0]  # A(0)->T(3), C(1)->G(2), G(2)->C(1), T(3)->A(0)

    # Precompute reverse complements for all possible stems of length k
    stem_rc_map = [0] * (4 ** k)
    for s in range(4 ** k):
        rc = 0
        t = s
        for _ in range(k):
            base_int = t % 4
            comp_base_int = comp_int[base_int]
            rc = rc * 4 + comp_base_int
            t = t // 4
        stem_rc_map[s] = rc

    # Convert the sequence to integer representation
    seq_int = [base_to_int.get(base, -1) for base in sequence]
    if -1 in seq_int:
        raise ValueError("Sequence contains invalid nucleotides. Only 'A', 'C', 'G', and 'T' are allowed.")

    # Precompute integer representations of all possible stems in the sequence
    max_pos = seq_len - k + 1
    stem_ints = [0] * max_pos
    s = 0
    for i in range(k):
        s = s * 4 + seq_int[i]
    stem_ints[0] = s

    for i in range(1, max_pos):
        s = ((s % (4 ** (k - 1))) * 4) + seq_int[i + k - 1]
        stem_ints[i] = s

    # Iterate through each possible stem start position
    for i in range(max_pos):
        # Calculate valid loop size positions
        start_j = i + k + min_loop
        end_j = min(i + k + max_loop + 1, max_pos)
        stem1_int = stem_ints[i]
        for j in range(start_j, end_j):
            stem2_int = stem_ints[j]
            stem2_rc_int = stem_rc_map[stem2_int]
            if stem1_int == stem2_rc_int:
                count += 1

    return count

# @njit
# def optimized_non_stupid_hairpin_counter(sequence, min_stem=3, min_loop=4, max_loop=9):
#     count = 0
#     seq_len = len(sequence)
#     # Map bases to integers
#     base_map = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
#     complement = np.array([3, 2, 1, 0], dtype=np.int8)  # A-T, C-G

#     # Convert sequence to array of integers
#     seq_array = np.empty(seq_len, dtype=np.int8)
#     for idx in range(seq_len):
#         seq_array[idx] = base_map.get(sequence[idx], -1)

#     # Precompute reverse complements for stems
#     num_stems = seq_len - min_stem + 1
#     stems = np.empty((num_stems, min_stem), dtype=np.int8)
#     reverse_complements = np.empty((num_stems, min_stem), dtype=np.int8)

#     for i in range(num_stems):
#         stem = seq_array[i:i + min_stem]
#         stems[i] = stem
#         reverse_complements[i] = complement[stem[::-1]]

#     # Iterate over possible stem starts
#     for i in range(seq_len - min_stem):
#         stem1 = stems[i]

#         min_j = i + min_stem + min_loop
#         max_j = min(i + min_stem + max_loop + 1, seq_len - min_stem + 1)

#         for j in range(min_j, max_j):
#             stem2 = stems[j]

#             # Compare stem1 with reverse complement of stem2
#             if np.array_equal(stem1, reverse_complements[j]):
#                 count += 1

#     return count

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
