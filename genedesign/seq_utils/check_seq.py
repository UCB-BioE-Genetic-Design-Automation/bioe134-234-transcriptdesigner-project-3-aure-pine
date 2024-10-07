from genedesign.checkers.codon_checker import CodonChecker
from genedesign.checkers.forbidden_sequence_checker import ForbiddenSequenceChecker
from genedesign.checkers.internal_promoter_checker import PromoterChecker
from genedesign.checkers.hairpin_checker import hairpin_checker
from genedesign.checkers.gc_content_checker import gc_content_checker

class CheckSequence:

    def __init__(self) -> None:
        self.codon_checker = None
        self.forbidden_checker = None
        self.promoter_checker = None
        
    def initiate(self) -> None:
        self.codon_checker = CodonChecker()
        self.forbidden_checker = ForbiddenSequenceChecker()
        self.promoter_checker = PromoterChecker()
        
        self.codon_checker.initiate()
        self.forbidden_checker.initiate()
        self.promoter_checker.initiate()

    def run(self, codons: list[str]) -> bool:
        results = []

        results.append(self.codon_checker.run(codons)[0])

        dna_seq = ''.join(codons)

        results.append(self.forbidden_checker.run(dna_seq)[0])
        results.append(self.promoter_checker.run(dna_seq)[0])
        results.append(hairpin_checker(dna_seq)[0])
        results.append(gc_content_checker(dna_seq))

        result = all(results)
        
        if not result:
            return False, results
        
        return result, None

