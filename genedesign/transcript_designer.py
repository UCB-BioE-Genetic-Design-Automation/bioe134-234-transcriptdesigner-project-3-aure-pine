from genedesign.rbs_chooser import RBSChooser
from genedesign.models.transcript import Transcript
from genedesign.checkers.codon_checker import CodonChecker
from genedesign.checkers.forbidden_sequence_checker import ForbiddenSequenceChecker
from genedesign.checkers.hairpin_checker import hairpin_checker
from genedesign.checkers.internal_promoter_checker import PromoterChecker
from seq_utils.sample_codon import SampleCodon
import numpy as np

class TranscriptDesigner:
    """
    Reverse translates a protein sequence into a DNA sequence and chooses an RBS using the highest CAI codon for each amino acid.
    """

    def __init__(self):
        self.aminoAcidToCodon = {}
        self.rbsChooser = None
        self.codon_checker = None
        self.forbidden_checker = None
        self.promoter_checker = None
        self.sample_codon = None


    def initiate(self) -> None:
        """
        Initializes the codon table and the RBS chooser.
        """
        ### FOR DEVELOPMENT
        self.rng = np.random.default_rng(seed=42)
        
        self.rbsChooser = RBSChooser()
        self.codon_checker = CodonChecker()
        self.forbidden_checker = ForbiddenSequenceChecker()
        self.promoter_checker = PromoterChecker()
        self.sample_codon = SampleCodon()

        self.rbsChooser.initiate()
        self.codon_checker.initiate()
        self.forbidden_checker.initiate()
        self.promoter_checker.initiate()
        self.sample_codon.initiate()

    def run(self, peptide: str, ignores: set) -> Transcript:
        

        
        #while index < len(peptide) -3:

        # Translate peptide to codons
    
        window_codons = NotImplemented
        selected_dnaseq = ''.join(window_codons)

        # Check sequence
        codon_bool = self.codon_checker(window_codons)[0]
        forbidden_bool = self.forbidden_checker(selected_dnaseq)[0]
        promoter_bool = self.promoter_checker(selected_dnaseq)[0]
        hairpin_bool = hairpin_checker(selected_dnaseq)[0]

        if not all([codon_bool, forbidden_bool, promoter_bool, hairpin_bool]):
            break

        index += 3
        
        window_codons = []

        for aa in window_peptides:
            codon = self.sample_codon.run(aa, self.rng)
            window_codons.append(codon)
        

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
