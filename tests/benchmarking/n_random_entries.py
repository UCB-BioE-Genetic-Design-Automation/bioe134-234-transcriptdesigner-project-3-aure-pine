import traceback
from genedesign.transcript_designer import TranscriptDesigner
import proteome_benchmarker as pb
import numpy as np
import time

def sample_entries(fasta_file, n_entries, rng):
    proteome = pb.parse_fasta(fasta_file)
    sample = rng.choice(list(proteome.items()), n_entries, replace=False, shuffle=False)
    return dict(sample)

def benchmark_proteome_sample(fasta_file, n_entries, rng):
    """
    Benchmarks the proteome using TranscriptDesigner.
    """
    designer = TranscriptDesigner()
    designer.initiate()

    proteome_sample = sample_entries(fasta_file, n_entries, rng)
    successful_results = []
    error_results = []

    for gene, protein in proteome_sample.items():
        try:
            print(f"Processing gene: {gene} with protein sequence: {protein[:30]}...")
            ignores = set()
            transcript = designer.run(protein, ignores)
            successful_results.append({
                'gene': gene,
                'protein': protein,
                'transcript': transcript
            })
        except Exception as e:
            error_results.append({
                'gene': gene,
                'protein': protein,
                'error': f"Error: {str(e)}\nTraceback: {traceback.format_exc()}"
            })
    
    return successful_results, error_results

def run_benchmark(fasta_file, n_entries, rng):
    """
    Runs the complete benchmark process: parsing, running TranscriptDesigner, validating, and generating reports.
    """
    start_time = time.time()
    
    # Benchmark the proteome
    parsing_start = time.time()
    successful_results, error_results = benchmark_proteome_sample(fasta_file, n_entries, rng)
    parsing_time = time.time() - parsing_start
    
    # Analyze and log errors
    errors_summary = pb.analyze_errors(error_results)
    
    # Validate the successful transcripts
    validation_start = time.time()
    validation_failures = pb.validate_transcripts(successful_results)
    execution_time = time.time() - validation_start

    # Write validation and error reports
    pb.write_validation_report(validation_failures)
    
    # Generate the summary report
    total_genes = len(successful_results) + len(error_results)
    pb.generate_summary(total_genes, parsing_time, execution_time, errors_summary, validation_failures)

if __name__ == "__main__":
    rng = np.random.default_rng(seed=42)
    entries_to_test = 5
    fasta_file = "tests/benchmarking/uniprotkb_proteome_UP000054015_2024_09_24.fasta"
    # print(sample_entries(fasta_file, entries_to_test, rng))
    run_benchmark(fasta_file, entries_to_test, rng) 