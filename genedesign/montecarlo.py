from genedesign.models.rbs_option import RBSOption
from genedesign.sliding_window_generator import sliding_window_generator
from genedesign.rbs_chooser import RBSChooser
from genedesign.seq_utils.sample_codon import SampleCodon
from genedesign.seq_utils.check_seq import CheckSequence
from genedesign.checkers.codon_checker import CodonChecker
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
        if not peptide:
            raise ValueError("Peptide needs to be a non-empty string.")

        codons = []
        full_peptide = peptide + '*' # Adding stop codon
        len_peptide = len(full_peptide)
        
        # Phase 1: Generate RBSOption
        first_6_aas = full_peptide[:6]
        first_6_codons = self.__montecarlo(first_6_aas, codons)
        gened_cds = ''.join(first_6_codons)
        selected_RBS = self.chooser.optimized_run(gened_cds, ignores)
        codons.extend(first_6_codons)

        # Phase 2:
        len_codons = len(codons)
        rest_peptide = full_peptide[len_codons:]

        for window in sliding_window_generator(rest_peptide, n_in_scope=self.n_codons_in_scope, n_ahead=self.n_ahead, step=self.step):
            window_codons = self.__find_codons(window, codons, selected_RBS, len_peptide)
            codons.extend(window_codons)

        return selected_RBS, codons
    
    def __find_codons(self, window, codons, selectedRBS, len_peptide):
            good_seq = False
            generated_codons = []
            attempts = 0
            max_attempts = 100  # Limit attempts to prevent infinite loop
            best_codon_score = 0
            best_generated_codons = []

            while not good_seq and attempts < max_attempts:
                # Generate
                # 3 (in scope) + n_ahead
                generated_codons = self.__montecarlo(window, codons)

                # Check
                good_seq, score = self.checker.run(generated_codons, codons, selectedRBS, len_peptide)

                # Compare
                if good_seq:
                    break
                elif score > best_codon_score:
                    best_codon_score = score
                    best_generated_codons = generated_codons

                attempts += 1

            if not good_seq:
                generated_codons = best_generated_codons
                # print(f'No valid window sequence found after {max_attempts}. Returning best codons.')

            return generated_codons[:self.n_codons_in_scope]

    def __montecarlo(self, window: str, last_n_codons: list[str]) -> list[str]:
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

        # Generate
        generated_codons = [self.sampler.run(amino_acid) for amino_acid in window]

        # Prepend the special first codon if applicable
        if special_first_codon:
            generated_codons = [special_first_codon] + generated_codons

        # Return the full window
        return generated_codons
    

    # def run1(self, peptide:str, ignores:set) -> tuple[RBSOption, list[str]]:
    #     """
    #     Phase 0:
    #     Generate first 6 codons for the RBSChooser, run RBSChooser and select a sequence
    #     # Choose an RBS
    #     selectedRBS = self.chooser.optimized_run(cds, ignores)

    #     PHASE 1:
    #     Generate up to N codons after those first 6
    #     Run the cds sequence tests, if those fail go back to the first 6 and generate again

    #     Phase 2:
    #     Window from N+1 to end, use the current window generation method to go to the end of the sequence with the tests

    #     Things to consider:
    #     How do I keep and store the N codons to check?
    #     How does the internal state
    #     Min length of input peptide: i know some are like 75 long
    #     """
    #     codons = []
    #     len_peptide = len(peptide)
    #     full_peptide = peptide + '*' # Adding stop codon

    #     # Phase 1:
    #     first_6_aas = full_peptide[:6]
    #     selectedRBS, first_6_codons = self.phase1(first_6_aas, ignores, codons, len_peptide)
    #     codons.extend(first_6_codons)

    #     # Phase 2:
    #     len(selectedRBS.utr) + 

    #     # Phase 3:
    #     rest_peptide = full_peptide[buffer_size:]
    #     for window in sliding_window_generator(rest_peptide, n_in_scope=self.n_codons_in_scope, n_ahead=self.n_ahead, step=self.step):
    #         window_codons = self.__montecarlo(window, codons, len_peptide, n_codons_in_scope=self.n_codons_in_scope)
    #         codons.extend(window_codons)
        
    #     return selectedRBS, codons

    # def oldrun(self, peptide:str, ignores:set) -> tuple[RBSOption, list[str]]:
    #     """
    #     Returns a tuple containing the generated RBSOption and the list of generated codons
    #     """
    #     # Find codon sequence
    #     codons = self.__find_codons(peptide)
    #     # Build the CDS from the codons
    #     cds = ''.join(codons)
    #     # Choose an RBS
    #     selectedRBS = self.chooser.optimized_run(cds, ignores)

    #     return selectedRBS, codons

# def __find_codons(self, peptide:str) -> list[str]:
#         codons = []
#         len_peptide = len(peptide)
#         full_peptide = peptide + '*' # Adding stop codon

#         for window in sliding_window_generator(full_peptide, n_in_scope=self.n_codons_in_scope, n_ahead=self.n_ahead, step=self.step):
#             # last_n_codons = codons[-self.n_behind:]
#             window_codons = self.__montecarlo(window, codons, len_peptide, n_codons_in_scope=self.n_codons_in_scope)
#             codons.extend(window_codons)

#         return codons
    

# # Hack solution
# # if len(last_n_codons) >= 31:
# #     good_seq, codon_diversity, rare_codon_count, cai_value = self.codon_checker.run(codons_to_check)
# # elif len(last_n_codons) < 31:
# #     good_seq = True


# # print(f"Good seq: {good_seq}, CAI: {cai_value}, Diversity: {codon_diversity}, Rare codons: {rare_codon_count}")



# if good_seq:
#     break

# current_metrics = (cai_value, codon_diversity, -rare_codon_count) # negative of rare codons to only use > in comparison

# if current_metrics > best_metrics: # comparison goes left to right, so left is most optimized side (cai, diversity, rare codons)
#     best_metrics = current_metrics
#     best_generated_codons = generated_codons
