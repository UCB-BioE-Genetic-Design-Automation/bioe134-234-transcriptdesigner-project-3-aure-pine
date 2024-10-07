from Bio import SeqIO
from collections import defaultdict

# Function to extract UTR, gene, and CDS information from the GenBank file
def extract_genes_info(genbank_file):
    gene_dict = defaultdict(dict)  # Dictionary to store gene info
    for record in SeqIO.parse(genbank_file, "genbank"):
        for feature in record.features:
            if feature.type == "gene":
                locus_tag = feature.qualifiers.get("locus_tag", [None])[0]
                gene_name = feature.qualifiers.get("gene", [None])[0]

                # CDS information
                cds_feature = None
                for cds in record.features:
                    if cds.type == "CDS" and cds.qualifiers.get("locus_tag") == [locus_tag]:
                        cds_feature = cds
                        break

                if cds_feature:
                    start, end = cds_feature.location.start, cds_feature.location.end
                    strand = cds_feature.location.strand
                    if strand == 1:  # Forward strand
                        utr_start = max(0, start - 50)
                        utr_seq = record.seq[utr_start:start]
                    else:  # Reverse strand, we need to reverse complement
                        utr_start = end
                        utr_seq = record.seq[utr_start:utr_start + 50].reverse_complement()

                    cds_seq = cds_feature.extract(record.seq)
                    # Save the gene information in the dictionary
                    gene_dict[locus_tag] = {
                        "gene": gene_name,
                        "UTR": utr_seq,
                        "CDS": cds_seq
                    }
    return gene_dict

def get_top_5_percent(file_path):
    """
    Extracts and returns the top 5% of tag:abundance pairs from a given text file.

    This function reads a text file containing tab-delimited lines of 'tag' and 'abundance'
    pairs. It strips the "511145." prefix from the tag and uses the remaining portion as the locus tag.
    It sorts the pairs by abundance in descending order and extracts the top 5% of entries.

    Args:
        file_path (str): The path to the input text file. The file should contain tab-delimited
                         lines where each line consists of a 'tag' (prefixed by "511145.") and
                         a corresponding 'abundance' value.

    Returns:
        dict: A dictionary where the keys are the 'locus_tag' values (strings, without the "511145." prefix)
              and the values are the corresponding 'abundance' (floats). This dictionary contains
              only the top 5% of the highest abundance entries.

    Raises:
        TypeError: If the file_path is not a string.
        FileNotFoundError: If the file at the specified path does not exist.
    """

    if not isinstance(file_path, str):
        raise TypeError("file_path must be a string representing the path to the input file.")

    data = []

    try:
        with open(file_path, 'r') as file:
            for line in file:
                # Skip comment lines and empty lines
                if line.startswith("#") or not line.strip():
                    continue

                # Split the line to extract the tag and abundance
                line_data = line.split()
                if len(line_data) != 2:
                    continue  # Skip invalid lines

                tag, abundance = line_data
                abundance = float(abundance)  # Convert abundance to float

                # Remove "511145." prefix from the tag
                if tag.startswith("511145."):
                    tag = tag.split("511145.")[1]

                data.append((tag, abundance))

    except FileNotFoundError:
        raise FileNotFoundError(f"The file '{file_path}' does not exist or cannot be found.")

    # Sort the data by abundance in descending order
    data.sort(key=lambda x: x[1], reverse=True)

    # Calculate the number of entries that make up the top 5%
    top_5_percent_count = int(len(data) * 0.05)

    # Get the top 5% of entries
    top_5_percent = data[:top_5_percent_count]

    return dict(top_5_percent)

def get_top_5_percent_utr_cds(locus_file_path, genbank_file_path):
    """
    Combines the top 5% of locus tags with their corresponding UTR, CDS, and gene name sequences.

    This function first extracts the top 5% of locus tags based on abundance from the text file.
    Then, it extracts the gene name, UTR (50 bp upstream), and CDS from the GenBank file for those
    top 5% locus tags.

    Args:
        locus_file_path (str): Path to the text file with locus tags and abundances.
        genbank_file_path (str): Path to the GenBank file containing sequence info.

    Returns:
        dict: A dictionary where keys are locus tags, and values are dictionaries containing:
              - 'gene': The gene name.
              - 'UTR': The upstream untranslated region (50 bp upstream of the CDS).
              - 'CDS': The coding sequence for that locus tag.
    """
    # Step 1: Get the top 5% locus tags from the text file
    top_5_percent_tags = get_top_5_percent(locus_file_path)

    # Step 2: Extract all gene info from the GenBank file
    gene_info = extract_genes_info(genbank_file_path)

    # Step 3: Filter the gene info for only the top 5% locus tags
    top_5_percent_info = {}
    for locus_tag in top_5_percent_tags:
        if locus_tag in gene_info:
            # Add gene, UTR, and CDS to the result
            top_5_percent_info[locus_tag] = {
                "gene": gene_info[locus_tag]["gene"],
                "UTR": gene_info[locus_tag]["UTR"],  # Using UTR (formerly "RBS")
                "CDS": gene_info[locus_tag]["CDS"]
            }

    return top_5_percent_info

# Example usage
locus_file_path = "genedesign/data/511145-WHOLE_ORGANISM-integrated.txt"  # Path to the text file
genbank_file_path = "genedesign/data/Ecoli_sequence.gb"  # Path to the GenBank file


top_5_percent_utr_cds = get_top_5_percent_utr_cds(locus_file_path, genbank_file_path)

# Print a few examples from the resulting dictionary
if top_5_percent_utr_cds:
    for locus_tag, sequences in list(top_5_percent_utr_cds.items())[:5]:  # Print first 5 entries
        print(f"Locus Tag: {locus_tag}")
        print(f"Gene: {sequences['gene']}")
        print(f"UTR: {sequences['UTR']}")
        print(f"CDS: {sequences['CDS']}")
        print("-" * 40)
else:
    print("No matching locus tags found in the top 5%.")