import subprocess
import sys
import os
from collections import defaultdict

#get the command-line arguments
bed_file = sys.argv[1]
fasta_file = sys.argv[2]
output_dir = sys.argv[3]

output_file = os.path.join(output_dir, "longest_transcript_and_features.tsv")
cmd_getfasta = f"bedtools getfasta -tab -bed {bed_file} -fi {fasta_file} -nameOnly -s > {output_file}"
subprocess.call(cmd_getfasta, shell=True)
metrics_file = os.path.join(output_dir, "nucleotide_metrics.tsv")

#calcualte nt proportion
def calculate_nucleotide_proportions(sequence):
    length = len(sequence)
    if length == 0:
        return {"A": 0, "C": 0, "T": 0, "G": 0}

    proportions = {
        "A": sequence.count("A") / length,
        "C": sequence.count("C") / length,
        "T": sequence.count("T") / length,
        "G": sequence.count("G") / length,
    }
    return proportions

sequences_by_header = defaultdict(str)

with open(output_file, "r") as f:
    for line in f:
        entry_name, sequence = line.strip().split("\t")
        sequences_by_header[entry_name] += sequence

with open(metrics_file, "w") as mf:

    mf.write("Feature\tA\tC\tT\tG\n")
    

    for header, combined_sequence in sequences_by_header.items():
        proportions = calculate_nucleotide_proportions(combined_sequence)
        mf.write(f"{header}\t{proportions['A']:.4f}\t{proportions['C']:.4f}\t{proportions['T']:.4f}\t{proportions['G']:.4f}\n")

