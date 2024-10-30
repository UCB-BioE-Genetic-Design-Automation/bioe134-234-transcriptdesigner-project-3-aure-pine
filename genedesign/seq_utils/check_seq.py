from genedesign.checkers.forbidden_sequence_checker import ForbiddenSequenceChecker
from genedesign.checkers.internal_promoter_checker import PromoterChecker
from genedesign.checkers.hairpin_checker import hairpin_checker

class CheckSequence:

    def __init__(self) -> None:
        self.codon_checker = None
        self.forbidden_checker = None
        self.promoter_checker = None
        
    def initiate(self) -> None:
        self.forbidden_checker = ForbiddenSequenceChecker()
        self.promoter_checker = PromoterChecker()
        
        self.forbidden_checker.initiate()
        self.promoter_checker.initiate()

    def run(self, codons: list[str]) -> bool:
        results = []
        dna_seq = ''.join(codons)

        results.append(self.forbidden_checker.run(dna_seq)[0])
        results.append(self.promoter_checker.run(dna_seq)[0])
        results.append(hairpin_checker(dna_seq)[0])
        """
        NEED TO ADD THE OTHER CHECKERS
        """

        result = all(results)
        
        if not result:
            return False, results
        
        return True, None