from numpy import inf
from genedesign.models.rbs_option import RBSOption
from genedesign.sliding_window_generator import sliding_window_generator
from genedesign.rbs_chooser import RBSChooser
from genedesign.seq_utils.sample_codon import SampleCodon
from genedesign.seq_utils.check_seq import CheckSequence
from genedesign.checkers.codon_checker import CodonChecker

# from genedesign.checkers.gc_content_checker import gc_content_checker
# from genedesign.checkers.codon_checker import CodonChecker
class MonteCarlo():
    def __init__(self):
        # Params
        self.n_codons_in_scope = 3
        self.n_behind = 3
        self.n_ahead = 6
        self.step = None

    def initiate(self):
        self.sampler = SampleCodon()
        self.chooser = RBSChooser()
        self.checker = CheckSequence()
        self.codon_checker = CodonChecker()

        self.sampler.initiate()
        self.chooser.initiate()
        self.checker.initiate()
        self.codon_checker.initiate()
    
    def run(self, peptide:str, ignores:set) -> tuple[RBSOption, list[str]]:
        """
        Returns a tuple containing the generated RBSOption and the list of generated codons
        """
        # Reset codon usages
        self.sampler.reset_codon_usages()
        # Find codon sequence
        codons = self.__find_codons(peptide)
        # Build the CDS from the codons
        cds = ''.join(codons)
        # Choose an RBS
        selectedRBS = self.chooser.optimized_run(cds, ignores)

        # good_seq = False
        # attempts = 0
        # max_attempts = 50  # Limit attempts to prevent infinite loop

        # while not good_seq and attempts < max_attempts:
        #     # Find codon sequence
        #     codons = self.__find_codons(peptide)
        #     # Build the CDS from the codons
        #     cds = ''.join(codons)
        #     # Choose an RBS
        #     selectedRBS = self.chooser.optimized_run(cds, ignores)

        #     full_seq = (selectedRBS.utr + cds).upper()
        #     good_seq = self.checker.run(full_seq)[0]
        #     attempts += 1

        return selectedRBS, codons
    
    def __find_codons(self, peptide:str) -> list[str]:
        codons = []
        full_peptide = peptide + '*' # Adding stop codon

        for window in sliding_window_generator(full_peptide, n_in_scope=self.n_codons_in_scope, n_ahead=self.n_ahead, step=self.step):
            # last_n_codons = codons[-self.n_behind:]
            window_codons = self.__montecarlo(window, codons, n_codons_in_scope=self.n_codons_in_scope)
            codons.extend(window_codons)

        return codons
    
    def __montecarlo(self, window: str, last_n_codons: list[str], n_codons_in_scope=3) -> list[str]:
        """
        Description:
        Generates a sequence of codons based on a window of amino acids and checks if the generated codons form a valid sequence using a sequence checker. 
        The function repeatedly generates codons until a valid sequence is found.
        """

        # Ensure the window is not empty
        if not window:
            raise ValueError("Window cannot be empty")
        
        # Handle special cases for the first codon if at the start of the sequence
        special_first_codon = None
        if not last_n_codons:
            first_amino_acid = window[0]
            special_codons = {'V': 'GTG', 'L': 'TTG', 'M': 'ATG'}
            special_first_codon = special_codons.get(first_amino_acid)

            # Exclude the first amino acid since it's handled
            window = window[1:]

        good_seq = False
        generated_codons = []
        attempts = 0
        max_attempts = 5000  # Limit attempts to prevent infinite loop

        # Initialize variables to store the best solution
        best_metrics = (0, 0, inf)
        best_generated_codons = None

        # Continue sampling codons and checking the sequence until a valid one is generated, or choose the one the minimizes 
        while not good_seq and attempts < max_attempts:

            generated_codons = [self.sampler.run(amino_acid) for amino_acid in window]
            # Prepend the special first codon if applicable
            if special_first_codon:
                generated_codons = [special_first_codon] + generated_codons
            
            codons_to_check = last_n_codons + generated_codons
            good_seq, codon_diversity, rare_codon_count, cai_value = self.codon_checker.run(codons_to_check)
            # print(f"Good seq: {good_seq}, CAI: {cai_value}, Diversity: {codon_diversity}, Rare codons: {rare_codon_count}")

            if good_seq:
                break

            current_metrics = (cai_value, codon_diversity, -rare_codon_count)

            if current_metrics > best_metrics: # comparison goes left to right, so left is most optimized side
                best_metrics = current_metrics
                best_generated_codons = generated_codons

            attempts += 1

        if not good_seq:
            generated_codons = best_generated_codons
            print(f"Unable to generate a valid sequence for a window after {max_attempts} attempts. Returning the best sequence found.")
        
        # Return only the first n codons in scope
        return generated_codons[:n_codons_in_scope]