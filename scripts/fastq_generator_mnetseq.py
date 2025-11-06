import pandas as pd
import gzip
import random
import argparse
import os
import subprocess
import shutil
import string
import ast
import numpy as np

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))

def run_cmd(cmd):
    print(f"[Running] {cmd}")
    subprocess.run(cmd, shell=True, check=True)


def make_random_suffix(length=5):
    return ''.join(random.choices(string.ascii_lowercase, k=length))


def extract_sequence_regions(df, read_length):
    df = df.copy()
    # Ensure converted_positions and incorporated_positions are lists of ints
    # Extract sequence after fragmentation
    # Extract read sequences

    df['read_upstream'] = df.apply(
        lambda row: row['parsed_sequence'][:read_length] if len(row['parsed_sequence']) >= read_length else row['parsed_sequence'],
        axis=1
    )
    df['read_downstream'] = df.apply(
        lambda row: row['parsed_sequence'][-read_length:] if len(row['parsed_sequence']) >= read_length else row['parsed_sequence'],
        axis=1
    )

    return df


def reverse_complement(seq):
    """Return the reverse complement of a DNA sequence."""
    complement = str.maketrans('ACGTacgt', 'TGCAtgca')
    return seq.translate(complement)[::-1]


def process_sequence(sequence, reverse_comp=False):
    """Process sequence by reversing and complementing if needed."""
    return reverse_complement(sequence) if reverse_comp else sequence


def convert_to_fastq(df, output_prefix, sequencing_type, strandedness, read_length, rates_df, ttseq=False):

    def get_sequences(row, which):
        # Always return a list of sequences
        seqs = row[which]
        if isinstance(seqs, list):
            return seqs
        elif isinstance(seqs, str):
            return [seqs]
        else:
            return []

    def get_absolute_coords(fragment_start_pos, fragment_end_pos):
        subset = rates_df[
            (rates_df['nucleotide_coord'] >= fragment_start_pos) &
            (rates_df['nucleotide_coord'] <= fragment_end_pos)
        ]
        if not subset.empty:
            chrom = subset['chromosome'].iloc[0]
            abs_start = subset['absolute_position'].iloc[0]
            abs_end = subset['absolute_position'].iloc[-1]
            strand = subset['strand'].iloc[0]
            return f"{chrom}:{abs_start}-{abs_end}:{strand}"
        else:
            return "NA"

    if sequencing_type == "SE":
        output_file = f"{output_prefix}.fastq.gz"
        with gzip.open(output_file, 'wt') as f:
            for _, row in df.iterrows():
                if strandedness == "rf":
                    seqs = get_sequences(row, 'read_downstream')
                    which_read = "downstream"
                    rev = True
                    abs_coords = get_absolute_coords(max(0, row['fragment_end_pos'] - read_length), row['fragment_end_pos'] - 1)
                elif strandedness == "fr":
                    seqs = get_sequences(row, 'read_upstream')
                    which_read = "upstream"
                    rev = False
                    abs_coords = get_absolute_coords(row['fragment_start_pos'], row['fragment_start_pos'] + read_length - 1)
                elif strandedness == "unstranded":
                    if random.choice([True, False]):
                        seqs = get_sequences(row, 'read_downstream')
                        which_read = "downstream"
                        rev = True
                        abs_coords = get_absolute_coords(max(0, row['fragment_end_pos'] - read_length), row['fragment_end_pos'] - 1)
                    else:
                        seqs = get_sequences(row, 'read_upstream')
                        which_read = "upstream"
                        rev = False
                        abs_coords = get_absolute_coords(row['fragment_start_pos'], row['fragment_start_pos'] + read_length - 1)

                for sequence in seqs:
                    seq = process_sequence(sequence, reverse_comp=rev)
                    random_suffix = make_random_suffix()
                    read_name = (
                        f"{row['molecule_id']}{random_suffix}"
                        f"_{abs_coords}"
                    )

                    quality_scores = "I" * len(seq)
                    f.write(f"{read_name}\n{seq}\n+\n{quality_scores}\n")

    elif sequencing_type == "PE":
        output_file_r1 = f"{output_prefix}_R1.fastq.gz"
        output_file_r2 = f"{output_prefix}_R2.fastq.gz"
        with gzip.open(output_file_r1, 'wt') as f1, gzip.open(output_file_r2, 'wt') as f2:
            for _, row in df.iterrows():
                upstream_seqs = get_sequences(row, 'read_upstream')
                downstream_seqs = get_sequences(row, 'read_downstream')
                n_pairs = min(len(upstream_seqs), len(downstream_seqs))

                for i in range(n_pairs):
                    random_suffix = make_random_suffix()  # single suffix per pair

                    if strandedness == "rf":
                        seq_r1 = process_sequence(downstream_seqs[i], reverse_comp=True)
                        seq_r2 = process_sequence(upstream_seqs[i], reverse_comp=False)

                        which_read_r1 = "downstream"
                        which_read_r2 = "upstream"

                        abs_coords_r1 = get_absolute_coords(max(0, row['fragment_end_pos'] - read_length), row['fragment_end_pos'] - 1)
                        abs_coords_r2 = get_absolute_coords(row['fragment_start_pos'], row['fragment_start_pos'] + read_length - 1)

                    elif strandedness == "fr":
                        seq_r1 = process_sequence(upstream_seqs[i], reverse_comp=False)
                        seq_r2 = process_sequence(downstream_seqs[i], reverse_comp=True)

                        which_read_r1 = "upstream"
                        which_read_r2 = "downstream"

                        abs_coords_r1 = get_absolute_coords(row['fragment_start_pos'], row['fragment_start_pos'] + read_length - 1)
                        abs_coords_r2 = get_absolute_coords(max(0, row['fragment_end_pos'] - read_length), row['fragment_end_pos'] - 1)

                    elif strandedness == "unstranded":
                        if random.choice([True, False]):
                            seq_r1 = process_sequence(downstream_seqs[i], reverse_comp=True)
                            seq_r2 = process_sequence(upstream_seqs[i], reverse_comp=False)

                            which_read_r1 = "downstream"
                            which_read_r2 = "upstream"

                            abs_coords_r1 = get_absolute_coords(max(0, row['fragment_end_pos'] - read_length), row['fragment_end_pos'] - 1)
                            abs_coords_r2 = get_absolute_coords(row['fragment_start_pos'], row['fragment_start_pos'] + read_length - 1)
                        else:
                            seq_r1 = process_sequence(upstream_seqs[i], reverse_comp=False)
                            seq_r2 = process_sequence(downstream_seqs[i], reverse_comp=True)

                            which_read_r1 = "upstream"
                            which_read_r2 = "downstream"

                            abs_coords_r1 = get_absolute_coords(row['fragment_start_pos'], row['fragment_start_pos'] + read_length - 1)
                            abs_coords_r2 = get_absolute_coords(max(0, row['fragment_end_pos'] - read_length), row['fragment_end_pos'] - 1)

                    random_suffix = make_random_suffix()
                    read_name_r1 = (
                        f"{row['molecule_id']}{random_suffix}_{abs_coords_r1}"
                    )
                    read_name_r2 = (
                        f"{row['molecule_id']}{random_suffix}_{abs_coords_r2}"
                    )
                    quality_scores_r1 = "I" * len(seq_r1)
                    quality_scores_r2 = "I" * len(seq_r2)
                    f1.write(f"{read_name_r1}\n{seq_r1}\n+\n{quality_scores_r1}\n")
                    f2.write(f"{read_name_r2}\n{seq_r2}\n+\n{quality_scores_r2}\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="convert dataframe to FASTQ file.")
    parser.add_argument("--input_df", help="full path to the input dataframe")
    parser.add_argument("--insert_size", type=str, default="200,300", help="Length of insert")
    parser.add_argument("--read_length", type=int, default=100, help="length of each read")
    parser.add_argument("--seq_type", choices=["PE", "SE"], type=str, default='SE', help="Type of sequencing: SE or PE")
    parser.add_argument("--s", choices=["rf", "fr", "unstranded"], type=str, default='rf', help="strandedness of the simulated library")
    parser.add_argument('--seq_depth', type=int, default=20000000, help='total library sequencing depth')
    parser.add_argument("--threads", type=int, default=1, help="Number of threads")
    parser.add_argument('--tpm_lower_limit', type=int, default=5, help='lower possible CPM per gene')
    parser.add_argument('--tpm_upper_limit', type=int, default=200, help='higher possible CPM per gene')
    parser.add_argument('--o', type=str, default='./', help='output path')
    parser.add_argument("--seed", type=int, help="random seed for reproducibility")

    args = parser.parse_args()

    if args.seed is not None:
        random.seed(args.seed)

    filename = os.path.splitext(os.path.basename(args.input_df))[0]
    base_filename = filename.split("_")[0]
    output_prefix = f"{args.o}/reads/{base_filename}"
    chopped_coordinates_file_path = f"{args.o}/mnetseqmRNAs/{base_filename}_mnetseq.tsv.gz"
    mRNA_coordinates_file_path = f"{args.o}/mRNA/{base_filename}.tsv.gz"

    try:
        df = pd.read_csv(chopped_coordinates_file_path, delimiter="\t")
        
        #Bring in original mRNA lengths to calculate tpms as for other methods
        df_mrna = pd.read_csv(mRNA_coordinates_file_path, delimiter="\t")
        df_mrna["sequence_length"] = df_mrna["full_molecule_sequence"].str.len()
        # mean of sequence_length column (divided by 1000)
        gene_length = df_mrna["sequence_length"].mean() / 1000

        # random uniform value between tpm_lower_limit and tpm_upper_limit
        gene_tpm = np.random.uniform(
            low=int(args.tpm_lower_limit),
            high=int(args.tpm_upper_limit)
            )
        # sequencing depth (divided by 1e6)
        seq_depth = args.seq_depth / 1e6
        # number of fragments for the gene
        reads_to_get = gene_length * gene_tpm * seq_depth

        if len(df) >= reads_to_get:
            df = df.sample(
                n=int(reads_to_get), replace=False, random_state=None
                )
        else:
            sampled_reads = df.sample(
                n=int(reads_to_get) - len(df), replace=True, random_state=None
                )
            # combine original and sampled
            df = pd.concat([df, sampled_reads], ignore_index=True)

        result_df = extract_sequence_regions(df, args.read_length)

        rates_for_gene = f"{args.o}/rate_per_gene/{base_filename}_RatesandTraversalTimes.gtf"
        rates_df = pd.read_csv(rates_for_gene, delimiter="\t")
        convert_to_fastq(result_df, output_prefix, sequencing_type=args.seq_type, strandedness=args.s,
                         read_length=args.read_length, rates_df=rates_df)
        #convert_to_fastq(result_df, output_prefix, sequencing_type=args.seq_type, strandedness=args.s, gene_id=base_filename)
        #try:
        #    os.remove(chopped_coordinates_file_path)
        #except OSError:
        #    pass
            
    except FileNotFoundError:
        pass
