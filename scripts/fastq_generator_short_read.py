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
    for col in ['converted_positions', 'incorporated_positions']:
        df[col] = df[col].apply(lambda x: ast.literal_eval(x) if isinstance(x, str) else x)

    # Explode read_coordinates
    df_exploded = df.copy()
    df_exploded['read_coordinates'] = df_exploded['read_coordinates'].str.split(',')
    df_exploded = df_exploded.explode('read_coordinates').reset_index()

    # Extract start/end coordinates
    coord_split = df_exploded['read_coordinates'].str.split('-', expand=True)
    df_exploded['read_start'] = coord_split[0].astype(int)
    df_exploded['read_end'] = coord_split[1].astype(int)

    # Extract read sequences
    df_exploded['read_upstream'] = df_exploded.apply(
        lambda row: df.loc[row['index'], 'full_molecule_sequence'][row['read_start']:row['read_start'] + read_length],
        axis=1
    )
    df_exploded['read_downstream'] = df_exploded.apply(
        lambda row: df.loc[row['index'], 'full_molecule_sequence'][max(0, row['read_end'] - read_length):row['read_end']],
        axis=1
    )

    # Count both converted and incorporated positions in upstream/downstream
    def count_and_list_subs(row):
        idx = row['index']
        read_start = row['read_start']
        read_end = row['read_end']

        # Regions
        upstream_start = read_start
        upstream_end = read_start + read_length
        downstream_start = max(0, read_end - read_length)
        downstream_end = read_end

        def extract_hits(pos_list):
            pos_list = [int(pos) for pos in pos_list]
            upstream_hits = [pos - upstream_start for pos in pos_list if upstream_start <= pos <= upstream_end]
            downstream_hits = [pos - downstream_start for pos in pos_list if downstream_start <= pos <= downstream_end]
            return upstream_hits, downstream_hits

        # Converted
        conv = df.loc[idx, 'converted_positions']
        conv_upstream_hits, conv_downstream_hits = extract_hits(conv)

        # Incorporated
        inc = df.loc[idx, 'incorporated_positions']
        inc_upstream_hits, inc_downstream_hits = extract_hits(inc)

        return pd.Series({
            'n_converted_in_upstream': len(conv_upstream_hits),
            'converted_positions_upstream': conv_upstream_hits,
            'n_converted_in_downstream': len(conv_downstream_hits),
            'converted_positions_downstream': conv_downstream_hits,
            'n_incorporated_in_upstream': len(inc_upstream_hits),
            'incorporated_positions_upstream': inc_upstream_hits,
            'n_incorporated_in_downstream': len(inc_downstream_hits),
            'incorporated_positions_downstream': inc_downstream_hits
        })

    df_exploded[
        ['n_converted_in_upstream', 'converted_positions_upstream',
         'n_converted_in_downstream', 'converted_positions_downstream',
         'n_incorporated_in_upstream', 'incorporated_positions_upstream',
         'n_incorporated_in_downstream', 'incorporated_positions_downstream']
    ] = df_exploded.apply(count_and_list_subs, axis=1)

    return df_exploded


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

    def get_absolute_coords(read_start, read_end):
        subset = rates_df[
            (rates_df['nucleotide_coord'] >= read_start) &
            (rates_df['nucleotide_coord'] <= read_end)
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
                if ttseq and not (
                    (row['n_incorporated_in_upstream'] > 0) or 
                    (row['n_incorporated_in_downstream'] > 0)
                ):
                    continue

                if strandedness == "rf":
                    seqs = get_sequences(row, 'read_downstream')
                    which_read = "downstream"
                    rev = True
                    abs_coords = get_absolute_coords(max(0, row['read_end'] - read_length), row['read_end'] - 1)
                elif strandedness == "fr":
                    seqs = get_sequences(row, 'read_upstream')
                    which_read = "upstream"
                    rev = False
                    abs_coords = get_absolute_coords(row['read_start'], row['read_start'] + read_length - 1)
                elif strandedness == "unstranded":
                    if random.choice([True, False]):
                        seqs = get_sequences(row, 'read_downstream')
                        which_read = "downstream"
                        rev = True
                        abs_coords = get_absolute_coords(max(0, row['read_end'] - read_length), row['read_end'] - 1)
                    else:
                        seqs = get_sequences(row, 'read_upstream')
                        which_read = "upstream"
                        rev = False
                        abs_coords = get_absolute_coords(row['read_start'], row['read_start'] + read_length - 1)

                for i, sequence in enumerate(seqs):
                    seq = process_sequence(sequence, reverse_comp=rev)
                    random_suffix = make_random_suffix()

                    if which_read == "upstream":
                        n_converted = row['n_converted_in_upstream']
                        converted_positions = row['converted_positions_upstream']
                        n_incorporated = row['n_incorporated_in_upstream']
                        incorporated_positions = row['incorporated_positions_upstream']
                    else:
                        n_converted = row['n_converted_in_downstream']
                        converted_positions = row['converted_positions_downstream']
                        n_incorporated = row['n_incorporated_in_downstream']
                        incorporated_positions = row['incorporated_positions_downstream']

                    conv_str = ','.join(str(pos + 1) for pos in converted_positions)
                    inc_str = ','.join(str(pos + 1) for pos in incorporated_positions)

                    read_name = (
                        f"{row['molecule_id']}{random_suffix}"
                        f"_{abs_coords}"
                        f"_ninc{n_incorporated}:{inc_str}_nsubs{n_converted}:{conv_str}"
                    )

                    quality_scores = "I" * len(seq)
                    f.write(f"{read_name}\n{seq}\n+\n{quality_scores}\n")

    elif sequencing_type == "PE":
        output_file_r1 = f"{output_prefix}_R1.fastq.gz"
        output_file_r2 = f"{output_prefix}_R2.fastq.gz"
        with gzip.open(output_file_r1, 'wt') as f1, gzip.open(output_file_r2, 'wt') as f2:
            for _, row in df.iterrows():
                if ttseq and not (
                    (row['n_incorporated_in_upstream'] > 0) or 
                    (row['n_incorporated_in_downstream'] > 0)
                ):
                    continue

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

                        abs_coords_r1 = get_absolute_coords(max(0, row['read_end'] - read_length), row['read_end'] - 1)
                        abs_coords_r2 = get_absolute_coords(row['read_start'], row['read_start'] + read_length - 1)

                    elif strandedness == "fr":
                        seq_r1 = process_sequence(upstream_seqs[i], reverse_comp=False)
                        seq_r2 = process_sequence(downstream_seqs[i], reverse_comp=True)

                        which_read_r1 = "upstream"
                        which_read_r2 = "downstream"

                        abs_coords_r1 = get_absolute_coords(row['read_start'], row['read_start'] + read_length - 1)
                        abs_coords_r2 = get_absolute_coords(max(0, row['read_end'] - read_length), row['read_end'] - 1)

                    elif strandedness == "unstranded":
                        if random.choice([True, False]):
                            seq_r1 = process_sequence(downstream_seqs[i], reverse_comp=True)
                            seq_r2 = process_sequence(upstream_seqs[i], reverse_comp=False)

                            which_read_r1 = "downstream"
                            which_read_r2 = "upstream"

                            abs_coords_r1 = get_absolute_coords(max(0, row['read_end'] - read_length), row['read_end'] - 1)
                            abs_coords_r2 = get_absolute_coords(row['read_start'], row['read_start'] + read_length - 1)
                        else:
                            seq_r1 = process_sequence(upstream_seqs[i], reverse_comp=False)
                            seq_r2 = process_sequence(downstream_seqs[i], reverse_comp=True)

                            which_read_r1 = "upstream"
                            which_read_r2 = "downstream"

                            abs_coords_r1 = get_absolute_coords(row['read_start'], row['read_start'] + read_length - 1)
                            abs_coords_r2 = get_absolute_coords(max(0, row['read_end'] - read_length), row['read_end'] - 1)

                    if which_read_r1 == "upstream":
                        n_conv_r1 = row['n_converted_in_upstream']
                        pos_conv_r1 = row['converted_positions_upstream']
                        n_inc_r1 = row['n_incorporated_in_upstream']
                        pos_inc_r1 = row['incorporated_positions_upstream']
                    else:
                        n_conv_r1 = row['n_converted_in_downstream']
                        pos_conv_r1 = row['converted_positions_downstream']
                        n_inc_r1 = row['n_incorporated_in_downstream']
                        pos_inc_r1 = row['incorporated_positions_downstream']

                    if which_read_r2 == "upstream":
                        n_conv_r2 = row['n_converted_in_upstream']
                        pos_conv_r2 = row['converted_positions_upstream']
                        n_inc_r2 = row['n_incorporated_in_upstream']
                        pos_inc_r2 = row['incorporated_positions_upstream']
                    else:
                        n_conv_r2 = row['n_converted_in_downstream']
                        pos_conv_r2 = row['converted_positions_downstream']
                        n_inc_r2 = row['n_incorporated_in_downstream']
                        pos_inc_r2 = row['incorporated_positions_downstream']

                    read_name_r1 = (
                        f"{row['molecule_id']}{random_suffix}_{abs_coords_r1}"
                        f"_ninc{n_inc_r1}:{','.join(str(pos + 1) for pos in pos_inc_r1)}"
                        f"_nsubs{n_conv_r1}:{','.join(str(pos + 1) for pos in pos_conv_r1)}"
                    )

                    read_name_r2 = (
                        f"{row['molecule_id']}{random_suffix}_{abs_coords_r2}"
                        f"_ninc{n_inc_r2}:{','.join(str(pos + 1) for pos in pos_inc_r2)}"
                        f"_nsubs{n_conv_r2}:{','.join(str(pos + 1) for pos in pos_conv_r2)}"
                    )

                    quality_scores_r1 = "I" * len(seq_r1)
                    quality_scores_r2 = "I" * len(seq_r2)

                    f1.write(f"{read_name_r1}\n{seq_r1}\n+\n{quality_scores_r1}\n")
                    f2.write(f"{read_name_r2}\n{seq_r2}\n+\n{quality_scores_r2}\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="convert dataframe to FASTQ file.")
    parser.add_argument("--input_df", help="full path to the input dataframe")
    parser.add_argument("--insert_size", type=str, default="200,300", help="Length of insert")
    parser.add_argument("--read_length", type=int,default=100,help="length of each read")
    parser.add_argument("--seq_type", choices=["PE", "SE"], type=str, default='SE', help="Type of sequencing: SE or PE")
    parser.add_argument("--s", choices=["rf", "fr","unstranded"], type=str, default='rf', help="strandedness of the simulated library")
    parser.add_argument('--seq_depth', type=int, default=20000000, help='total library sequencing depth')
    parser.add_argument("--threads", type=int, default=1, help="Number of threads")
    parser.add_argument('--tpm_lower_limit', type=int, default=5, help='lower possible CPM per gene')
    parser.add_argument('--tpm_upper_limit', type=int, default=200, help='higher possible CPM per gene')
    parser.add_argument("--fragments", action="store_true", help="export a ground truth for fragmentation & size selection")
    parser.add_argument('--o', type=str, default='./', help='output path')
    parser.add_argument("--seed", type=int, help="random seed for reproducibility")

    args = parser.parse_args()

    if args.seed is not None:
        random.seed(args.seed)

    cmd_chop = [
            f"Rscript {os.path.join(SCRIPT_DIR, 'short_read_chopper.R')}",
            f"--tsv {args.input_df}",
            f"--insert_size {args.insert_size}",
            f"--read_length {args.read_length}",
            f"--threads {args.threads}",
            f"--seq_depth {args.seq_depth}",
            f"--tpm_lower_limit {args.tpm_lower_limit}",
            f"--tpm_upper_limit {args.tpm_upper_limit}",
            f"-o {args.o}"
        ]
    if args.seed:
        cmd_chop += ['--seed', args.seed]
    if args.fragments:
        cmd_chop += ['--fragments with_ground_truth']

    run_cmd(" ".join(cmd_chop))

    filename = os.path.splitext(os.path.basename(args.input_df))[0]
    base_filename = filename.split(".")[0]
    output_prefix = f"{args.o}/reads/{base_filename}"
    chopped_coordinates_file_path = f"{args.o}/temp/mRNAs_with_fragments/{base_filename}_fragments.tsv"
    rates_for_gene = f"{args.o}/rate_per_gene/{base_filename}_RatesandTraversalTimes.gtf"

    try:
        df = pd.read_csv(chopped_coordinates_file_path, delimiter="\t")
        #Only keep rows (mRNAs that have at least one incorporated nucleotide analog, not doing this on the chopper because it would remove all background)
        df = df[df['incorporated_positions'].map(bool)]

        for col in ['converted_positions', 'incorporated_positions']:
            df[col] = df[col].apply(lambda x: ast.literal_eval(x) if isinstance(x, str) else (x if isinstance(x, list) else []))

        rates_df = pd.read_csv(rates_for_gene, delimiter="\t")
        result_df = extract_sequence_regions(df, args.read_length)

        #Add background if any
        path_to_BGmRNAs=os.path.join(os.path.dirname(args.o), "mRNA", base_filename + "_background.tsv.gz")
        if os.path.exists(path_to_BGmRNAs) and os.path.getsize(path_to_BGmRNAs) > 0:
            cmd_chop = [
            f"Rscript {os.path.join(SCRIPT_DIR, 'short_read_chopper.R')}",
            f"--tsv {path_to_BGmRNAs}",
            f"--insert_size {args.insert_size}",
            f"--read_length {args.read_length}",
            f"--threads {args.threads}",
            f"-o {args.o}"
            ]
            if args.seed:
                cmd_chop += ['--seed', args.seed]

            run_cmd(" ".join(cmd_chop))
            #Same as the pulldown samples but skipping the enrichment since these molecules are not labeled
            chopped_coordinates_BG_file_path = f"{args.o}/temp/mRNAs_with_fragments/{base_filename}_background_fragments.tsv"
            df_bg = pd.read_csv(chopped_coordinates_BG_file_path, delimiter="\t")
            df_bg['molecule_id'] = df_bg['molecule_id'].astype(str) + '_BG_'#Add to the molecule ID that this is a background molecule
            for col in ['converted_positions', 'incorporated_positions']:
                df_bg[col] = df_bg[col].apply(lambda x: ast.literal_eval(x) if isinstance(x, str) else (x if isinstance(x, list) else []))



            result_df_bg = extract_sequence_regions(df_bg, args.read_length)
            result_df = pd.concat([result_df, result_df_bg], ignore_index=True)

        


        #Bring in original mRNA lengths to calculate tpms as for other methods
        df_mrna = pd.read_csv(args.input_df, delimiter="\t")
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

        if len(result_df) >= reads_to_get:
            result_df = result_df.sample(
                n=int(reads_to_get), replace=False, random_state=None
                )
        else:
            sampled_reads = result_df.sample(
                n=int(reads_to_get) - len(result_df), replace=True, random_state=None
                )
            # combine original and sampled
            result_df = pd.concat([result_df, sampled_reads], ignore_index=True)



        convert_to_fastq(result_df, output_prefix, sequencing_type=args.seq_type, strandedness=args.s,
                         read_length=args.read_length, rates_df=rates_df)
        try:
            os.remove(chopped_coordinates_file_path)
        except OSError:
            pass
    except FileNotFoundError:
        pass





