import argparse
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from BCBio import GFF
import re

# Define a function to count motif occurrences
# Important, this function is able to detect motifs that wrap around a newline in standard fasta format. Therefore, if you manually search the motifs with grep or a text editor, you may find fewer than this function does - but the function works correctly.
def count_motif(sequence, motif):
    motif_regex = motif.replace('N', '.')  # Replace N with . for any nucleotide match
    return len(re.findall(motif_regex, str(sequence)))

# Define a function to calculate expected motif frequency
def calculate_expected_frequency(motif, gc_content):
    # Convert the percentage GC content to a probability (between 0 and 1)
    gc_frequency = gc_content / 100
    # Calculate the AT frequency as the remainder of the GC frequency subtracted from 1
    at_frequency = 1 - gc_frequency
    # Ensure that the motif is in uppercase for consistent comparison
    motif = motif.upper()
    # Initialize the probability of the entire motif occurring to 1 (100%)
    motif_probability = 1
    # Iterate over each base in the motif
    for base in motif:
        # If the base is A or T, multiply the existing motif_probability by the probability of A or T
        if base in ['A', 'T']:
            motif_probability *= at_frequency / 2
        # If the base is G or C, multiply the existing motif_probability by the probability of G or C
        elif base in ['G', 'C']:
            motif_probability *= gc_frequency / 2
        # If the base is 'N' (which stands for any nucleotide), it does not change the probability
        elif base == 'N':
            motif_probability *= 1
    # Multiply by 1000 for 1 kb and by 2 for the two strands
    expected_frequency_per_kb = motif_probability * 1000 * 2
    return expected_frequency_per_kb

# Set up the argument parser
parser = argparse.ArgumentParser(description='Count motifs in gene upstream regions and write upstream sequences to FASTA.')
parser.add_argument('-f', '--fasta', required=True, help='Path to the FASTA genome file.')
parser.add_argument('-a', '--gff', required=True, help='Path to the GFF annotation file.')
parser.add_argument('-g', '--genes', required=True, help='Path to the file containing gene IDs.')
parser.add_argument('-l', '--upstream_length', type=int, required=True, help='Length of the upstream region to search for motifs in basepairs.')
parser.add_argument('-p', '--prefix', type=str, required=True, help='Prefix for output filenames')
parser.add_argument('-m', '--motifs', help='Path to the file containing motifs.')
parser.add_argument('-cm', '--count-motifs', action='store_true', help='Flag to trigger counting of motifs.')
parser.add_argument('-fq', '--frequency', action='store_true', help='Flag to trigger calculation of motif frequencies (expected and observed).')
parser.add_argument('-gc', '--gc_content', type=float, required=False, help='GC content of the genome (percentage).')

# Parse the arguments
args = parser.parse_args()

# Load the FASTA file as a dictionary
fasta_sequences = SeqIO.to_dict(SeqIO.parse(args.fasta, "fasta"))

# Read the gene IDs from file
with open(args.genes, 'r') as f:
    gene_ids = [line.strip() for line in f]

# Read the motifs from file if counting is enabled
motifs = {}
if args.count_motifs:
    if args.motifs:
        with open(args.motifs, 'r') as m:
            for line in m:
                parts = line.strip().split("\t")
                if len(parts) == 2:
                    motifs[parts[0]] = parts[1]
    else:
        raise ValueError("Motif file must be provided if --cm is TRUE (-m or --motifs).")

    # Prepare the output table header
    header = ['Gene_ID'] + list(motifs.keys())

    # Prepare the frequencies file
    if args.frequency:
        frequencies_output_file_name = f"{args.prefix}_frequencies.txt"
        frequencies_output_file = open(frequencies_output_file_name, 'w')
        frequencies_output_file.write("\t".join(['Motif', 'Expected', 'Observed']) + "\n")

# Check if motif counting is enabled and open the output file for writing motif counts
if args.count_motifs:
    counts_output_file_name = f"{args.prefix}_motif_counts.{args.upstream_length // 1000}kb_upstream.txt"
    counts_output_file = open(counts_output_file_name, 'w')
    counts_output_file.write("\t".join(header) + "\n")

if args.frequency:
    # Initialize the dictionary to count motifs across all genes
    motif_counts_across_genes = {motif_name: 0 for motif_name in motifs.keys()}

# Open the output file for writing upstream sequences
fasta_output_file_name = f"{args.prefix}_upstream_sequences_{args.upstream_length // 1000}kb.fasta"
with open(fasta_output_file_name, 'w') as fasta_output_file:
    # Parse the GFF file
    with open(args.gff, "r") as gff_file:
        for record in GFF.parse(gff_file, base_dict=fasta_sequences):
            for feature in record.features:
                if feature.type == "gene":
                    gene_id = feature.qualifiers.get('ID', [None])[0]
                    for sub_feature in feature.sub_features:
                        if sub_feature.type == "mRNA":
                            mrna_id = sub_feature.qualifiers.get('ID', [None])[0]
                            if mrna_id in gene_ids:
                                gene_start = feature.location.start.position
                                gene_end = feature.location.end.position
                                strand = feature.strand

                                if strand == 1:  # '+' strand
                                    start = max(0, gene_start - args.upstream_length)
                                    end = gene_start
                                else:  # '-' strand
                                    start = gene_end
                                    end = min(len(record.seq), gene_end + args.upstream_length)

                                # Extract the sequence
                                upstream_seq = record.seq[start:end]
                                if strand == -1:
                                    upstream_seq = upstream_seq.reverse_complement()

                                # Create a SeqRecord
                                upstream_record = SeqRecord(upstream_seq, id=mrna_id, description='')

                                # Write the upstream sequence to the FASTA file
                                SeqIO.write(upstream_record, fasta_output_file, "fasta")

                                # Count motifs if enabled
                                if args.count_motifs:
                                    motif_counts = [mrna_id]
                                    for motif_name, motif_seq in motifs.items():
                                        motif_count = count_motif(upstream_seq, motif_seq)
                                        motif_counts.append(str(motif_count))
                                        # If user enables calculation of motif frequency across genes
                                        if args.frequency:
                                            # Add the motif count to the tally across all genes
                                            motif_counts_across_genes[motif_name] += motif_count
                                    # Write the counts for this gene to the counts output file
                                    counts_output_file.write("\t".join(motif_counts) + "\n")



# Close the counts file if motif counting was enabled
if args.count_motifs:
    counts_output_file.close()
    print("Motif counting complete. Counts saved to:", counts_output_file_name)

# Calculate and write expected and observed frequencies to frequencies file
if args.frequency:
    total_upstream_length = args.upstream_length * len(gene_ids)
    total_length_kb = total_upstream_length / 1000
    for motif_name, motif_seq in motifs.items():
        print(motif_name+": "+motif_seq)
        print("Counts: "+str(motif_counts_across_genes[motif_name]))
        total_motif_count = motif_counts_across_genes[motif_name]
        # frequency per kb!
        expected_frequency = calculate_expected_frequency(motif_seq, args.gc_content)
        observed_frequency = total_motif_count / total_length_kb
        # Calculate the percentage difference
        percentage_difference = ((observed_frequency - expected_frequency) / expected_frequency) * 100
        # Write the results to the output file
        frequencies_output_file.write("\t".join([
            motif_name,
            "{:.6f}".format(expected_frequency),
            "{:.6f}".format(observed_frequency),
            "{:.2f}%".format(percentage_difference)
        ]) + "\n")

    frequencies_output_file.close()

    print("Motif frequencies calculated. Frequencies saved to:", frequencies_output_file_name)

print("Upstream sequences saved to:", fasta_output_file_name)
