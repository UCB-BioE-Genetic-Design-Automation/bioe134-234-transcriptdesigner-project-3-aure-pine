from genedesign.models.rbs_option import RBSOption
from genedesign.seq_utils.calc_edit_distance import calculate_edit_distance
from genedesign.seq_utils.hairpin_counter import non_stupid_hairpin_counter, optimized_non_stupid_hairpin_counter
from genedesign.seq_utils.Translate import Translate
from genedesign.seq_utils.get_top_5_percent_utr_cds import get_top_5_percent_utr_cds
import time
import Levenshtein

class RBSChooser:
    """
    A simple RBS selection algorithm that chooses an RBS from a list of options, excluding any RBS in the ignore set.
    """

    def __init__(self):
        self.translator = None
        self.rbsOptions = []

    def initiate(self) -> None:
        """
        Initialization method for RBSChooser.
        """

        locus_file_path = "genedesign/data/511145-WHOLE_ORGANISM-integrated.txt"  # Path to the text file
        genbank_file_path = "genedesign/data/Ecoli_sequence.gb"  # Path to the GenBank file

        self.translator = Translate()
        self.translator.initiate()
        self.rbs_options: list[RBSOption] = []

        top_5_percent_utr_cds = get_top_5_percent_utr_cds(locus_file_path, genbank_file_path)

        for locus_tag, data in top_5_percent_utr_cds.items():
            gene_name = data['gene']
            utr = data['UTR']
            cds = data['CDS']

            # Calculate the first six amino acids
            first_six_aas = self.translator.run(cds[:18])  # First six amino acids come from the first 18 nucleotides

            # Create an RBSOption and add it to the list
            rbs_option = RBSOption(utr=utr, cds=cds, gene_name=gene_name, first_six_aas=first_six_aas)
            self.rbs_options.append(rbs_option)


    def run(self, cds: str, ignores: set[RBSOption]) -> RBSOption:
        """
        Executes the RBS selection process for the given CDS.

        Parameters:
        - cds (str): The coding sequence to pair with an RBS. Must contain only valid DNA nucleotides (A, T, G, C).
        - ignores (Set[RBSOption]): A set of RBSOption instances to ignore during selection.

        Returns:
        - RBSOption: The selected RBSOption that best pairs with the given CDS.

        Raises:
        - ValueError: If no valid RBS options remain after filtering.
        """
        # Step 1: Exclude RBSOptions that are in the ignores set
        available_rbs_options = [rbs_option for rbs_option in self.rbs_options if rbs_option not in ignores]

        # Check if any valid RBSOptions remain after exclusion
        if not available_rbs_options:
            raise ValueError("No valid RBS options available after filtering.")

        # Initialize variables to keep track of the best option and the best score
        best_rbs_option = None
        best_score = float('inf')  # Lower score is better
        input_peptide = self.translator.run(cds[:18])  # Get the first six amino acids from the input CDS

        # Define weights for each criterion (adjust these as needed)
        hairpin_weight = 0.5  # Weight for secondary structure
        peptide_weight = 0.5  # Weight for peptide similarity

        total_hairpin_time = 0
        total_similarity_time = 0
        # Step 2 and 3: Iterate through available RBSOptions
        for rbs_option in available_rbs_options:
            hc_start = time.time()
            # Step 2: Calculate the hairpin score using hairpin_counter
            hairpin_score = optimized_non_stupid_hairpin_counter(rbs_option.utr + cds, min_stem=4, min_loop=3, max_loop=8)
            total_hairpin_time += time.time() - hc_start

            psc_start = time.time()
            # Step 3: Calculate the peptide similarity score using calculate_edit_distance
            peptide_similarity_score = calculate_edit_distance(input_peptide, rbs_option.first_six_aas)
            total_similarity_time += time.time() - psc_start

            # Combine the scores using a weighted sum
            final_score = (float(hairpin_score) * hairpin_weight) + (float(peptide_similarity_score) * peptide_weight)

            # Update the best RBSOption if the current one has a lower score
            if final_score < best_score:
                best_score = final_score
                best_rbs_option = rbs_option

        print(f"    Total hairpin_counter time: {total_hairpin_time}")
        print(f'    Total similarity time:      {total_similarity_time}')

        # Return the RBSOption with the lowest combined score
        return best_rbs_option
    
    def optimized_run(self, cds: str, ignores: set[RBSOption]) -> RBSOption:
        """
        Executes the RBS selection process for the given CDS.

        Parameters:
        - cds (str): The coding sequence to pair with an RBS.
        - ignores (Set[RBSOption]): A set of RBSOption instances to ignore during selection.

        Returns:
        - RBSOption: The selected RBSOption that best pairs with the given CDS.
        """
        # Exclude ignored RBSOptions
        available_rbs_options = [
            rbs_option for rbs_option in self.rbs_options if rbs_option not in ignores
        ]

        if not available_rbs_options:
            raise ValueError("No valid RBS options available after filtering.")

        input_peptide = self.translator.run(cds[:18])  # First six amino acids

        hairpin_weight = 0.5
        peptide_weight = 0.5

        best_score = float('inf')
        best_rbs_option = None

        for rbs_option in available_rbs_options:
            # Hairpin score
            hairpin_score = optimized_non_stupid_hairpin_counter(
                rbs_option.utr + cds, min_stem=4, min_loop=3, max_loop=8
            )

            # Peptide similarity score
            peptide_similarity_score = Levenshtein.distance(
                input_peptide, rbs_option.first_six_aas
            )

            # Final weighted score
            final_score = (hairpin_score * hairpin_weight) + (peptide_similarity_score * peptide_weight)

            # Update the best option
            if final_score < best_score:
                best_score = final_score
                best_rbs_option = rbs_option

        return best_rbs_option

if __name__ == "__main__":
    # Example usage of RBSChooser
    cds1 = "ATGGTAAGAAAACAGTTGCAGAGAGTTGAATT"
    cds2 = "ATGGCTAGCAAATACGATTTTACAATATAA"

    # Initialize the chooser
    chooser = RBSChooser()
    chooser.initiate()

    # Choose RBS with no ignores
    ignores = set()
    selected1 = chooser.run(cds2, ignores)
    
    # Add the first selection to the ignore list
    ignores.add(selected1)
    
    # Choose another RBS option after ignoring the first
    selected2 = chooser.run(cds2, ignores)

    # Print the selected RBS options
    print("Selected1:", selected1)
    print("Selected2:", selected2)
