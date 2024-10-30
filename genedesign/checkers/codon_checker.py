import sys
import csv
from collections import Counter  # Import Counter for counting codons

class CodonChecker:
    """
    Description: 
    This class checks codon usage for a given CDS and calculates metrics like:
    1. Codon Diversity: Fraction of unique codons.
    2. Rare Codon Count: Number of rare codons in the sequence.
    3. Codon Adaptation Index (CAI): Based on codon usage frequencies.

    Input (run method):
    cds (List[str]): A list of codons representing the coding sequence (e.g., ['ATG', 'TAA', 'CGT']).

    Output:
    Tuple[bool, float, int, float]: A tuple containing:
        - codons_above_board (bool): True if the CDS passes the thresholds, False otherwise.
        - codon_diversity (float): Fraction of unique codons.
        - rare_codon_count (int): Number of rare codons.
        - cai_value (float): CAI value.
    """

    codon_frequencies: dict[str, float]
    rare_codons: list[str]
    rare_codon_threshold: float

    def initiate(self) -> None:
        """
        Loads codon usage data from a file and sets up the codon frequencies and rare codons.
        """
        codon_usage_file = 'genedesign/data/codon_usage.txt'
        self.codon_frequencies = {}
        self.rare_codons = []
        self.rare_codon_threshold = 0.1  # Threshold for rare codon frequency

        with open(codon_usage_file, 'r') as f:
            reader = csv.reader(f, delimiter='\t')
            for row in reader:
                if len(row) < 4:
                    continue  # Skip invalid rows
                codon = row[0].strip()
                usage_freq = float(row[2].strip())
                self.codon_frequencies[codon] = usage_freq
                
                # Identify rare codons
                if usage_freq < self.rare_codon_threshold:
                    self.rare_codons.append(codon)

    def run(self, cds: list[str]) -> tuple[bool, float, int, float]:
        """
        Calculates codon diversity, rare codon count, and Codon Adaptation Index (CAI) for the provided CDS.
        Returns a boolean indicating whether the codons pass specified thresholds.

        :param cds: List of codons representing the CDS.
        :return: Tuple containing a boolean, codon diversity, rare codon count, and CAI score.
        """
        if not cds:
            return False, 0.0, 0, 0.0  # Return false for empty CDS

        # Calculate codon diversity
        codon_counts = Counter(cds)
        total_codons = len(cds)
        codon_diversity = len(codon_counts) / total_codons if total_codons > 0 else 0.0

        # Count rare codons
        rare_codon_count = sum(codon_counts[codon] for codon in self.rare_codons if codon in cds)

        # Calculate CAI (Codon Adaptation Index) as the geometric mean of codon frequencies
        cai_numerators = [self.codon_frequencies.get(codon, 0.01) for codon in cds]  # Use 0.01 for unknown codons
        cai_product = 1
        for freq in cai_numerators:
            cai_product *= freq
        
        cai_value = cai_product ** (1 / len(cai_numerators)) if cai_numerators else 0.0

        # Apply thresholds to determine if the codons are above board
        diversity_threshold = 0.5
        rare_codon_limit = 3
        cai_threshold = 0.2

        codons_above_board = (codon_diversity >= diversity_threshold and
                              rare_codon_count <= rare_codon_limit and
                              cai_value >= cai_threshold)

        return codons_above_board, codon_diversity, rare_codon_count, cai_value

class CodonChecker2:
    """
    Description:
    This class checks codon usage for a given CDS and calculates metrics like:
    1. Codon Diversity: Fraction of unique codons.
    2. Rare Codon Count: Number of rare codons in the sequence.
    3. Codon Adaptation Index (CAI): Based on codon usage frequencies.
    """

    def initiate(self) -> None:
        """
        Loads codon usage data from a file and sets up the codon frequencies and rare codons.
        """
        codon_usage_file = 'genedesign/data/codon_usage.txt'
        self.codon_frequencies = {}       # {amino_acid: {codon: frequency}}
        self.max_codon_frequencies = {}   # {amino_acid: max_frequency}
        self.codon_to_amino_acid = {}     # {codon: amino_acid}
        self.rare_codons = []
        self.rare_codon_threshold = 0.1  # Threshold for rare codon frequency

        with open(codon_usage_file, 'r') as f:
            reader = csv.reader(f, delimiter='\t')
            for row in reader:
                if len(row) < 4:
                    continue  # Skip invalid rows
                codon = row[0].strip()
                amino_acid = row[1].strip()
                usage_freq = float(row[2].strip())
                # prob = float(row[3].strip())  # If there's a probability column

                # Initialize dictionaries
                if amino_acid not in self.codon_frequencies:
                    self.codon_frequencies[amino_acid] = {}
                self.codon_frequencies[amino_acid][codon] = usage_freq
                self.codon_to_amino_acid[codon] = amino_acid

                # Update the maximum frequency for the amino acid
                if amino_acid not in self.max_codon_frequencies or usage_freq > self.max_codon_frequencies[amino_acid]:
                    self.max_codon_frequencies[amino_acid] = usage_freq

                # Identify rare codons
                if usage_freq < self.rare_codon_threshold:
                    self.rare_codons.append(codon)

    def run(self, cds: list[str]) -> tuple[bool, float, int, float]:
        """
        Calculates codon diversity, rare codon count, and Codon Adaptation Index (CAI) for the provided CDS.
        Returns a boolean indicating whether the codons pass specified thresholds.
        """
        if not cds:
            return False, 0.0, 0, 0.0  # Return false for empty CDS

        # Calculate codon diversity
        codon_counts = Counter(cds)
        total_codons = len(cds)
        codon_diversity = len(codon_counts) / total_codons if total_codons > 0 else 0.0

        # Count rare codons
        rare_codon_count = sum(codon_counts[codon] for codon in self.rare_codons if codon in cds)

        # Calculate CAI (Codon Adaptation Index)
        w_values = []
        for codon in cds:
            amino_acid = self.codon_to_amino_acid.get(codon)
            if amino_acid and amino_acid in self.max_codon_frequencies:
                max_freq = self.max_codon_frequencies[amino_acid]
                codon_freq = self.codon_frequencies[amino_acid].get(codon, 0.01)
                w_i = codon_freq / max_freq
            else:
                w_i = 0.01  # Assign a small value if amino acid or max_freq is not found
            w_values.append(w_i)

        # Calculate the geometric mean of w_values
        cai_product = 1
        for w in w_values:
            cai_product *= w

        cai_value = cai_product ** (1 / len(w_values)) if w_values else 0.0

        # Apply thresholds to determine if the codons are above board
        diversity_threshold = 0.5
        rare_codon_limit = 3
        cai_threshold = 0.2

        codons_above_board = (codon_diversity >= diversity_threshold and
                              rare_codon_count <= rare_codon_limit and
                              cai_value >= cai_threshold)

        return codons_above_board, codon_diversity, rare_codon_count, cai_value

if __name__ == "__main__":
    """
    Main method for running the CodonChecker on a hardcoded CDS.
    This method initializes CodonChecker, runs it on a hardcoded CDS, and prints
    the results including a boolean indicating if the codons are acceptable,
    codon diversity, rare codon count, and CAI for the sequence.
    """
    # Initialize the codon checker
    codon_checker = CodonChecker()
    codon_checker.initiate()

    # Hardcoded example CDS
    cds = ['ATG', 'CAA', 'GGG', 'TAA']  # High CAI example
    # cds = ['AGG', 'AGA', 'AGG', 'AGA']  # Very low CAI example
    
    # Run CodonChecker
    codons_above_board, codon_diversity, rare_codon_count, cai_value = codon_checker.run(cds)

    # Print results
    print(f"CDS: {cds}")
    print(f"Codons Above Board: {codons_above_board}")
    print(f"Codon Diversity: {codon_diversity}")
    print(f"Rare Codon Count: {rare_codon_count}")
    print(f"Codon Adaptation Index (CAI): {cai_value}")
