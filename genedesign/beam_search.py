import heapq
from genedesign.models.rbs_option import RBSOption
from genedesign.rbs_chooser import RBSChooser
class BeamSearch():
    def __init__(self):
        ## heuristic func hyperparameters
        self.cai_weight = None
        self.hairpin_weight = None
        self.forbidden_seq_weight = None
        self.internal_promoter_weight = None
        self.rnase_weight = None

    def initiate(self):
        self.chooser = RBSChooser()
        self.chooser.initiate()

    def run(self, peptide:str, ignores:set) -> tuple[RBSOption, list[str]]:
        # Initialize the beam with sequences for the first amino acid
        # Iterate over the rest of the amino acid sequence
        #   Keep the top k sequences, using min priority queue? (min score is better, less bad stuff. Have to invert CAI tho)
        """
        Performs beam search to generate codon sequences from an amino acid sequence,
        with constraints checking and pruning of sequences that violate constraints.

        Parameters:
        - amino_acid_seq: list of amino acids (strings)
        - codon_table: dict mapping amino acids to lists of possible codons
        - beam_width: int, the number of sequences to keep at each step
        - codon_usage_table: dict mapping codons to their relative adaptiveness (w values)
        - min_cai_threshold: float, the minimum acceptable CAI score
        - forbidden_sequences: list of DNA sequences that are prohibited (e.g., RNAse binding sites)

        Returns:
        - List of codon sequences (each sequence is a list of codon strings)
        """
        # Constraint: The first amino acid-codon pair must be one of the specified pairs
        start_codons = {'V': 'GTG', 'L': 'TTG', 'M': 'ATG'}
        
        # Check if the first amino acid matches one of the required amino acids
        first_amino_acid = peptide[0]
        if first_amino_acid not in start_codons:
            raise ValueError(f"The first amino acid must be one of {list(start_codons.keys())}, but got '{first_amino_acid}'.")

        # Initialize the beam with the specified amino acid-codon pair
        initial_codon = start_codons[first_amino_acid]
        seq = [initial_codon]

        # Generate RBSOption
        rbs = self.__find_rbs(seq)

        score = self.calculate_score(seq, rbs)
        beam = [(score, seq, rbs)]

        for aa in peptide[1:]:
            candidates = []

            # Process current beam
            while beam:
                current_score, seq, utr_sequence = heapq.heappop(beam)
                possible_codons = codon_table.get(aa, [])
                if not possible_codons:
                    raise ValueError(f"No codons found for amino acid '{aa}'.")

                for codon in possible_codons:
                    new_seq = seq + [codon]

                    # Generate RBSOption
                    rbs = self.__find_rbs(new_seq)

                    # Check constraints
                    if self.violates_constraints(new_seq, forbidden_sequences, new_utr_sequence):
                        continue  # Prune this sequence

                    # Compute new score
                    new_score = compute_score(
                        new_seq, codon_usage_table, new_utr_sequence, mfe_weight, cai_weight, mfe_threshold
                    )

                    # Add to candidates
                    heapq.heappush(candidates, (new_score, new_seq, new_utr_sequence))

                    # Keep the beam size within the beam_width
                    if len(candidates) > beam_width:
                        heapq.heappop(candidates)  # Remove the worst candidate

            # Update the beam with the new candidates
            beam = candidates

            if not beam:
                raise Exception(
                    "No sequences left in the beam. Consider increasing the beam width or adjusting weights."
                )

        # After processing all amino acids, the best sequence is at the top of the heap
        best_score, best_seq, best_utr_sequence = heapq.heappop(beam)

        return best_utr_sequence, best_seq

    def __find_rbs(self, codons:list[str]) -> RBSOption:
        coding_sequence = ''.join(codons)
        return self.chooser.optimized_run(coding_sequence)

    def calculate_score(self, codon_sequence):
        # Compute the heuristic score for the sequence
        codon_score = NotImplemented # if using min, have to invert, or have it be inverted already here
        hairpin_score = NotImplemented
        fseq_score = NotImplemented
        ip_score = NotImplemented
        rnase_score = NotImplemented

        # might have to normalize here so scores can be compared easier
        # as sequences get longer, hairpin counts, # of internal promoters etc increase

        """
        length of longest rare codon stretch or num of rare codons?
        number of RNAase E sites

        coding region:
        - Codon usage (CAI)
        - internal starts
        - internal terminators
        - RNAse E sites

        first third OR -4 -> +37:
        - secondary structure

        5' UTR: (maybe don't need this one with RBSChooser)
        Shine-dalgarno
        secondary structure
        spacing
        """





        score = self.cai_weight * codon_score + self.hairpin_weight * hairpin_score + self.forbidden_seq_weight * fseq_score + self.internal_promoter_weight * ip_score + self.rnase_weight * rnase_score
        return score
    
    def violates_constraints(seq):
        # Check for forbidden sequences, RNAse sites, internal promoters, and other hard constraints
        # Return True if constraints are violated, False otherwise

        pass