import csv
import time
import traceback
from genedesign.checkers.codon_checker import CodonChecker
from genedesign.checkers.forbidden_sequence_checker import ForbiddenSequenceChecker
from genedesign.checkers.hairpin_checker import hairpin_checker
from genedesign.checkers.internal_promoter_checker import PromoterChecker
from genedesign.rbs_chooser import RBSChooser
from genedesign.seq_utils.Translate import Translate
from tests.benchmarking.proteome_benchmarker import generate_summary, validate_transcripts

result_file_path = '/Users/ian/1 - Projects/HW/bioe134/bioe134-234-transcriptdesigner-project-3-aure-pine/tests/benchmarking/genome_benchmark_results/'

def parse_fasta_gene_sequences(file_path):
    """
    Parses a FASTA file and returns a list of dictionaries with gene names and DNA sequences.

    Args:
        file_path (str): Path to the FASTA file.

    Returns:
        list: A list of dictionaries, each with 'gene' and 'transcript' keys.
    """
    gene_transcripts = []
    with open(file_path, 'r') as f:
        current_gene = None
        current_sequence = []
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                # If we have a current gene, save its sequence before moving to the next
                if current_gene:
                    gene_transcripts.append({
                        'gene': current_gene,
                        'transcript': ''.join(current_sequence)
                    })
                # Reset for the next gene
                current_sequence = []
                # Extract gene name from the header line using string methods
                start = line.find('[gene=')
                if start != -1:
                    start += len('[gene=')
                    end = line.find(']', start)
                    if end != -1:
                        current_gene = line[start:end]
                    else:
                        current_gene = None  # Handle the case where ']' is missing
                else:
                    current_gene = None  # Handle the case where '[gene=' is missing
            else:
                current_sequence.append(line)
        # Save the last gene sequence after the loop ends
        if current_gene:
            gene_transcripts.append({
                'gene': current_gene,
                'transcript': ''.join(current_sequence)
            })
    return gene_transcripts

# Example usage:
# gene_list = parse_fasta_gene_sequences('path_to_your_file.txt')
# print(gene_list)

def generate_genes(fasta_file):
    genes = parse_fasta_gene_sequences(fasta_file)

    rbs_chooser = RBSChooser()
    rbs_chooser.initiate()

    successful_results = []
    error_results = []

    for gene, sequence in genes.items():
        try:
            print(f"Processing gene: {gene} with protein sequence: {sequence[:30]}...")
            ignores = set()
            rbs = rbs_chooser.run(sequence, ignores)
            successful_results.append({
                'gene': gene,
                'sequence': sequence,
                'rbs': rbs
            })
        except Exception as e:
            error_results.append({
                'gene': gene,
                'sequence': sequence,
                'error': f"Error: {str(e)}\nTraceback: {traceback.format_exc()}"
            })
    
    return genes, error_results

def analyze_errors(error_results):
    """
    Write the error analysis to a text file.
    """
    error_summary = {}
    with open(f'{result_file_path}error_summary.txt', 'w') as f:
        for error in error_results:
            error_message = error['error'].split("\n")[0]
            error_summary[error_message] = error_summary.get(error_message, 0) + 1
            f.write(f"Gene: {error['gene']}\n{error['error']}\n\n")
    
    return error_summary

def validate_sequences(genes):
    """
    Validate the successful sequences using various checkers, now including CodonChecker.
    """
    forbidden_checker = ForbiddenSequenceChecker()
    forbidden_checker.initiate()
    promoter_checker = PromoterChecker()
    promoter_checker.initiate()
    translator = Translate()
    translator.initiate()
    codon_checker = CodonChecker()  # Initialize CodonChecker
    codon_checker.initiate()  # Load the codon usage data

    validation_failures = []
    for gene in genes:
        print(f'Validating {gene}')
        cds = gene['sequence']
        try:
            # Check if CDS length is a multiple of 3
            if len(cds) % 3 != 0:
                raise ValueError("CDS length is not a multiple of 3.")
            
            # Ensure CDS starts with valid start codon and ends with stop codon
            if not (cds.startswith(("ATG", "GTG", "TTG")) and cds.endswith(("TAA", "TGA", "TAG"))):
                raise ValueError("CDS does not start with a valid start codon or end with a valid stop codon.")
            
        except ValueError as e:
            validation_failures.append({
                'gene': gene,
                'cds': cds,
                'site': f"Translation or completeness error: {str(e)}"
            })
            continue

        # Validate against hairpins, forbidden sequences, and internal promoters
        transcript_dna = result['transcript'].rbs.utr.upper() + cds
        passed_hairpin, hairpin_string = hairpin_checker(transcript_dna)
        if not passed_hairpin:
            formatted_hairpin = hairpin_string.replace('\n', ' ').replace('"', "'")
            validation_failures.append({
                'gene': result['gene'],
                'protein': result['protein'],
                'cds': transcript_dna,
                'site': f"Hairpin detected: {formatted_hairpin}"
            })

        passed_forbidden, forbidden_site = forbidden_checker.run(transcript_dna)
        if not passed_forbidden:
            validation_failures.append({
                'gene': result['gene'],
                'protein': result['protein'],
                'cds': transcript_dna,
                'site': f"Forbidden sequence: {forbidden_site}"
            })

        passed_promoter, found_promoter = promoter_checker.run(transcript_dna)
        if not passed_promoter:
            validation_failures.append({
                'gene': result['gene'],
                'protein': result['protein'],
                'cds': transcript_dna,
                'site': f"Constitutive promoter detected: {found_promoter}" if found_promoter else "Constitutive promoter detected"
            })

        codons_above_board, codon_diversity, rare_codon_count, cai_value = codon_checker.run(result['transcript'].codons)
        if not codons_above_board:
            validation_failures.append({
                'gene': result['gene'],
                'protein': result['protein'],
                'cds': cds,
                'site': f"Codon usage check failed: Diversity={codon_diversity}, Rare Codons={rare_codon_count}, CAI={cai_value}"
            })
    
    return validation_failures

def write_validation_report(validation_failures):
    """
    Writes validation results to a TSV file.
    """
    with open(f'{result_file_path}validation_failures.tsv', 'w', newline='') as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerow(['gene', 'protein', 'cds', 'site'])
        for failure in validation_failures:
            writer.writerow([failure['gene'], failure['protein'], failure['cds'], failure['site']])

def run_benchmark(dna_fasta_file):
    """
    Runs the complete benchmark process: parsing, running TranscriptDesigner, validating, and generating reports.
    """
    parsing_start = time.time()
    sequences = generate_genes(dna_fasta_file)
    parsing_time = time.time() - parsing_start

    validation_start = time.time()
    validation_failures = validate_transcripts(sequences)
    execution_time = time.time() - validation_start

    # Write validation and error reports
    write_validation_report(validation_failures)
    
    # Generate the summary report
    total_genes = len(sequences)
    generate_summary(total_genes, parsing_time, execution_time, validation_failures)

if __name__ == "__main__":
    fasta_file = "genedesign/data/U00096_full_sequence.txt"
    run_benchmark(fasta_file)
