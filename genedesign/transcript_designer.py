from genedesign.models.transcript import Transcript

## SEARCH ALGORITHMS
from genedesign.montecarlo import MonteCarlo
from genedesign.beam_search import BeamSearch
class TranscriptDesigner:

    def __init__(self):
        self.search_algorithm = None

    def initiate(self) -> None:
        # Monte Carlo Search
        self.search_algorithm = MonteCarlo()
        self.search_algorithm.initiate()

        # Beam Search
        # self.search_algorithm = BeamSearch()
        # self.search_algorithm.initiate()

        # ML method??

    def run(self, peptide: str, ignores: set) -> Transcript:
        selectedRBS, codons = self.search_algorithm.run(peptide, ignores)

        # Return the Transcript object
        return Transcript(selectedRBS, peptide, codons)

if __name__ == "__main__":
    # Example usage of TranscriptDesigner
    GFP = "MSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTLVTTFSYGVQCFSRYPDHMKQHDFFKSAMPEGYVQERTIFFKDDGNYKTRAEVKFEGDTLVNRIELKGIDFKEDGNILGHKLEYNYNSHNVYIMADKQKNGIKVNFKIRHNIEDGSVQLADHYQQNTPIGDGPVLLPDNHYLSTQSALSKDPNEKRDHMVLLEFVTAAGITHGMDELYK"
    peptide = "MYPFIRTARMTV"
    
    designer = TranscriptDesigner()
    designer.initiate()

    ignores = set()
    transcript = designer.run(peptide, ignores)
    
    # Print out the transcript information
    print(transcript)
