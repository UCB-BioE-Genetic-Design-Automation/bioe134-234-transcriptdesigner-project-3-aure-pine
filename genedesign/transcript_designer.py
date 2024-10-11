from genedesign.rbs_chooser import RBSChooser
from genedesign.models.transcript import Transcript
from genedesign.seq_utils.check_seq import CheckSequence
from seq_utils.sample_codon import SampleCodon
from sliding_window_generator import sliding_window_generator
import numpy as np

## SEARCH ALGORITHMS
from genedesign.montecarlo import montecarlo

class TranscriptDesigner:
    """
    Reverse translates a protein sequence into a DNA sequence and chooses an RBS using the highest CAI codon for each amino acid.
    """

    def __init__(self):
        self.rbsChooser = None
        self.search_selection_algo = None
        self.seq_checker = None
        self.codon_sampler = None

        self.search_selection_algo = None

    def initiate(self) -> None:
        """
        Initializes the codon table and the RBS chooser.
        """
        ### FOR DEVELOPMENT
        self.rng = np.random.default_rng(seed=42)

        self.rbsChooser = RBSChooser()
        self.seq_checker = CheckSequence()
        self.codon_sampler = SampleCodon()

        self.rbsChooser.initiate()
        self.seq_checker.initiate()
        self.codon_sampler.initiate()

        ## SEARCH ALGORITHMS ##
        # Montecarlo method
        self.search_selection_algo = montecarlo

        # MCTS + window
        # self.search_selection_algo = MCTS()
        # self.search_selection_algo.initiate()

        # ML algo

    def run(self, peptide: str, ignores: set) -> Transcript:

        codons = []
        
        # Parameters
        n_codons_in_scope = 3
        n_behind = 3
        n_ahead = 6
        step = None

        for window in sliding_window_generator(peptide, n_in_scope=n_codons_in_scope, n_ahead=n_ahead, step=step):
            last_n_codons = codons[-n_behind:]
            window_codons = self.search_selection_algo.run(window, last_n_codons, n_codons_in_scope=n_codons_in_scope, codon_sampler=self.codon_sampler, seq_checker=self.seq_checker)
            codons.append(window_codons)

        # Build the CDS from the codons
        cds = ''.join(codons)

        # Choose an RBS
        selectedRBS = self.rbsChooser.run(cds, ignores)

        # Return the Transcript object
        return Transcript(selectedRBS, peptide, codons)

if __name__ == "__main__":
    # Example usage of TranscriptDesigner
    peptide = "MYPFIRTARMTV"
    
    designer = TranscriptDesigner()
    designer.initiate()

    ignores = set()
    transcript = designer.run(peptide, ignores)
    
    # Print out the transcript information
    print(transcript)
