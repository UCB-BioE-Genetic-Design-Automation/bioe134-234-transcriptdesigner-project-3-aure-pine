from collections import Counter
from genedesign.checkers.forbidden_sequence_checker import ForbiddenSequenceChecker
from genedesign.checkers.internal_promoter_checker import PromoterChecker
from genedesign.checkers.codon_checker import CodonChecker
from genedesign.checkers.hairpin_checker import hairpin_checker
from genedesign.models.rbs_option import RBSOption
from genedesign.checkers.gc_content_checker import gc_checker

class CheckSequence:

    def __init__(self) -> None:
        self.codon_checker = None
        self.forbidden_checker = None
        self.promoter_checker = None

        self.last_checked_positions = None
        self.window_sizes = None
        self.gc_count = None
        
    def initiate(self) -> None:
        self.forbidden_checker = ForbiddenSequenceChecker()
        self.promoter_checker = PromoterChecker()
        self.codon_checker = CodonChecker()
        
        self.forbidden_checker.initiate()
        self.promoter_checker.initiate()
        self.codon_checker.initiate()

        self.gc_count = 0

        # Initialize last checked positions for each checker
        self.last_checked_positions = {
            'forbidden': 0,
            'promoter': 0,
            'hairpin': 0,
            'gc_content': 0
        }
        
        # Define window sizes for each checker
        self.window_sizes = {
            'forbidden': 8,  
            'promoter': 29,    
            'hairpin': 50,                                    
            'gc_content': 50                                  
        }

    def run(self, generated_codons: list[str], codons: list[str], rbs: RBSOption, len_peptide) -> tuple[bool, float]:
        results = []
        
        all_codons = codons + generated_codons
        good_codons, cai = self.check_codons(all_codons, len_peptide)
        results.append(good_codons)

        #only want 25 bp into the rbs utr
        # 50 - (25 + 18, already generated) = 7, so window size should be 3? This is perfect
        
        dna_seq = ''.join(all_codons)
        full_seq = self.combine_sequences(rbs.utr, dna_seq)
        seq_length = len(full_seq)

        checker_results = []

        # ---- Forbidden Sequence Checker ----
        # Get last checked position and calculate new base pairs
        last_pos = self.last_checked_positions['forbidden']
        num_new_bp = seq_length - last_pos
        max_forbidden_seq_length = self.window_sizes['forbidden']

        if num_new_bp > 0:
            start_pos = max(0, last_pos - (max_forbidden_seq_length - 1))
            forbidden_seq_to_check = full_seq[start_pos:seq_length]
            result = self.forbidden_checker.run(forbidden_seq_to_check)[0]
            checker_results.append(result)
            self.last_checked_positions['forbidden'] = seq_length
        else:
            pass # can add different handling here

        # ---- GC Content Checker ----
        last_pos_gc = self.last_checked_positions['gc_content']
        num_new_bp_gc = seq_length - last_pos_gc
        gc_min_length = self.window_sizes['gc_content']

        if seq_length >= gc_min_length and num_new_bp_gc > 0:
            start_pos_gc = max(0, last_pos - (max_forbidden_seq_length - 1))
            sequence_to_check_gc = full_seq[start_pos_gc:seq_length]
            gc_result = gc_checker(sequence_to_check_gc)[0]
            checker_results.append(gc_result)
            self.last_checked_positions['gc_content'] = seq_length
        else:
            pass


        # ---- Other Checkers ----
        # Checkers with windows
        for checker_name in ['promoter', 'hairpin']:
            last_pos = self.last_checked_positions[checker_name]
            window_size = self.window_sizes[checker_name]
            
            # Calculate the number of new base pairs since last checked
            num_new_bp = seq_length - last_pos

            if num_new_bp >= window_size:
                # Adjust start position to include overlaps
                start_pos = max(0, last_pos - (window_size - 1))
                end_pos = seq_length
                
                # Extract the sequence to check
                sequence_to_check = full_seq[start_pos:end_pos]

                # Run the checker
                if checker_name == 'promoter':
                    result = self.promoter_checker.run(sequence_to_check)[0]
                elif checker_name == 'hairpin':
                    result = hairpin_checker(sequence_to_check)[0]

                # Append the result
                checker_results.append(result)

                # Update the last checked position
                self.last_checked_positions[checker_name] = seq_length
            else:
                # Not enough new base pairs; skip this checker
                # Optionally, you can assume the previous result remains valid
                pass  # Or handle as needed

        # Append checker results to overall results
        results.extend(checker_results)

        # Compute overall result and score
        num_true = sum(results)
        result = all(results)
        normalized_num_true = num_true / len(results)

        # Higher is better
        score = normalized_num_true + cai
        
        if not result:
            return False, score
        
        return True, score
    
    def combine_sequences(self, utr, cds):
        max_window_size = 50  # Desired total length of the output string
        max_chars_utr = 25  # Maximum characters to take from utr

        len_utr = len(utr)
        len_cds = len(cds)

        if len_cds >= max_window_size:
            # Case 2: cds is 50 characters or more
            result = cds[-max_window_size:]  # Take the last 50 characters from cds
        else:
            # Case 1: cds is less than 50 characters
            n = len_cds  # Number of characters to take from cds
            m = max_window_size - n  # Number of characters to take from utr

            # Adjust m if utr is shorter than m, at most max_chars_utr (25)
            m = min(m, len_utr, max_chars_utr)
            n = max_window_size - m  # Recalculate n in case m was adjusted

            # Extract the last m characters from utr (up to 25)
            last_m_chars_utr = utr[-m:] if m > 0 else ''
            # Extract the first n characters from cds
            first_n_chars_cds = cds[:n] if n > 0 else ''

            result = last_m_chars_utr + first_n_chars_cds

        return result
    
    def check_codons(self, codons, len_peptide):
        """
        Implemented logic here so codon_checker was not edited. I don't know how you are planning to run the tests and the checkers.
        """
        # diversity_weight = 1
        # cai_weight = 1
        # rare_weight = 1

        diversity_threshold = 0.5
        global_rare_codon_limit = 3
        cai_threshold = 0.2

        _, codon_diversity, rare_codon_count, cai_value = self.codon_checker.run(codons)

        num_codons = len(codons)

        # Codon diversity
        if num_codons < 62:
            num_diff_codons = int(codon_diversity * 62)
            actual_codon_diversity = num_diff_codons / num_codons
        else:
            actual_codon_diversity = codon_diversity
        
        # Finding the rare codon limit by location in the peptide
        # Divides the protein into 3 sequences, where += 1 rare codon can be used.
        rare_codon_limit = self.rare_codon_limit(num_codons, len_peptide, global_rare_codon_limit)

        good_seq = (actual_codon_diversity >= diversity_threshold and
                    rare_codon_count <= rare_codon_limit and
                    cai_value >= cai_threshold)

        return good_seq, cai_value


    def rare_codon_limit(self, current_length, peptide_length, global_rare_codon_limit):
        section_length = peptide_length // global_rare_codon_limit
        
        # Determine which section we're in based on current length
        section_index = min(current_length // section_length, global_rare_codon_limit - 1)

        # Calculate rare_codon_limit directly
        if section_index < global_rare_codon_limit - 1:
            rare_codon_limit = section_index + 1
        else:
            rare_codon_limit = global_rare_codon_limit

        return rare_codon_limit