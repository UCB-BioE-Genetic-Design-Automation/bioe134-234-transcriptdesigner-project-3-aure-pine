import numpy as np
import csv

class SampleCodon:
    """
    Description:
    A class that reads codon usage data from a file, stores it in a dictionary structure,
    and allows sampling codons based on Codon Adaptation Index (CAI) probabilities for a given amino acid.

    Input:
    - The codon usage data file (tab-separated) containing codon, amino acid, and CAI probabilities.
    - Amino acid input for sampling a codon.

    Output:
    - Codon sampled based on the CAI probabilities for the given amino acid.
    """
    
    def __init__(self) -> None:
        self.codon_probabilities = None
        self.rng = None
        self.amino_acids = None

    def initiate(self) -> None:
        """
        Reads codon usage data from a file, populates the dictionary structure,
        and converts the codon lists and probabilities to numpy arrays for efficient sampling.

        The dictionary structure:
        - Keys: Amino acid single-letter codes.
        - Values: A tuple of two numpy arrays:
            1. Array of codons for the corresponding amino acid.
            2. Array of CAI probabilities for each codon.

        Example:
        'A': (['GCT', 'GCC', 'GCA', 'GCG'], [0.4, 0.3, 0.2, 0.1])
        """
        self.rng = np.random.default_rng()

        self.amino_acids = np.array([
                'A', 'R', 'N', 'D', 'C', 
                'Q', 'E', 'G', 'H', 'I', 
                'L', 'K', 'M', 'F', 'P', 
                'S', 'T', 'W', 'Y', 'V',
                '*'
            ])

        # Initialize dictionary with empty lists for codons and probabilities
        self.codon_probabilities = {
            'A': ([], []), 'C': ([], []), 'D': ([], []), 'E': ([], []),
            'F': ([], []), 'G': ([], []), 'H': ([], []), 'I': ([], []),
            'K': ([], []), 'L': ([], []), 'N': ([], []), 'P': ([], []),
            'Q': ([], []), 'R': ([], []), 'S': ([], []), 'T': ([], []),
            'V': ([], []), 'W': ([], []), 'Y': ([], []), 
            'M': ([], []),  # Start codon
            '*': ([], [])  # Stop codons
        }

        # Add values to the dictionary from the file
        codon_usage_file = 'genedesign/data/codon_usage.txt'

        with open(codon_usage_file, 'r') as f:
            reader = csv.reader(f, delimiter='\t')
            for row in reader:
                if len(row) < 4:
                    continue  # Skip invalid rows
                codon = row[0].strip()        # Codon
                amino_acid = row[1].strip()   # Amino acid
                prob = float(row[2].strip())  # CAI value (probability)
                
                # Add the codon and its probability to the corresponding amino acid
                codons, probabilities = self.codon_probabilities[amino_acid]
                codons.append(codon)
                probabilities.append(prob)

        # Convert all lists into np.array for speed improvements
        for amino_acid in self.codon_probabilities:
            codons, probabilities = self.codon_probabilities[amino_acid]
            self.codon_probabilities[amino_acid] = (
                np.array(codons),           # Convert codons list to numpy array
                np.array(probabilities)     # Convert probabilities list to numpy array
            )

    def run(self, amino_acid:str) -> str:
        """
        Samples a codon for a given amino acid based on the stored CAI probabilities.

        Args:
        amino_acid (str): The amino acid for which a codon should be sampled.
                          The input should be a single-letter amino acid code (e.g., 'A', 'M', '*').

        Returns:
        str: The sampled codon based on the CAI probabilities for the given amino acid.
        """
        if not amino_acid in self.amino_acids:
            raise ValueError(f"Invalid amino acid: {amino_acid}.")
        codons, probs = self.codon_probabilities[amino_acid]
        return str(self.rng.choice(codons, p=probs))

    # For debugging
    def datastructure(self):
        return self.codon_probabilities
    
    def get_codons(self, amino_acid) -> list:
        return self.codon_probabilities[amino_acid][0]
    
    def get_usages(self, amino_acid) -> list:
        return self.codon_probabilities[amino_acid][1]
    
    def get_data(self, amino_acid) -> tuple[list, list]:
        return self.codon_probabilities[amino_acid]

if __name__ == "__main__":
    aa1 = "I"
    
    sampler = SampleCodon()
    sampler.initiate()

    ignores = set()
    codon1 = sampler.run(aa1)
    codon2 = sampler.run(aa1)
    
    # Print out the selected codons
    print(codon1), print(codon2)
    
    