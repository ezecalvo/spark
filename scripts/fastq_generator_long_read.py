import pandas as pd
import random
import string
import argparse
import gzip
import os

def generate_random_string(length=8):
    """Generate a random string of given length."""
    return ''.join(random.choices(string.ascii_letters + string.digits, k=length))

def reverse_complement(sequence):
    """Return the reverse complement of a DNA sequence."""
    complement = str.maketrans('ACGTacgt', 'TGCAtgca')
    return sequence.translate(complement)[::-1]

def process_sequence(sequence, sequencing_type):
    """Process sequence based on sequencing type."""
    if sequencing_type == "RNA":
        return sequence[::-1]  # Reverse the sequence
    elif sequencing_type == "cDNA":
        return sequence if random.choice([True, False]) else reverse_complement(sequence)
    return sequence

def id_generator(size=12, chars=string.ascii_uppercase + string.digits):
    return "@" + ''.join(random.choice(chars) for _ in range(size - 1))   

def convert_to_fastq(df, output_file, sequencing_type,base_filename):
    with gzip.open(output_file, 'wt') as f:
        for _, row in df.iterrows():
            converted_positions = row['converted_positions']
            nconvreads = len(converted_positions)
            converted_positions = row['converted_positions']
            processed_positions = converted_positions.strip('[]')
            processed_positions = processed_positions.replace(', ', '-')
            molecule_id=id_generator()

            read_name = f"{molecule_id}_{base_filename}"#_nconvbases:{nconvreads}_posconvbases:{processed_positions}"
            sequence = process_sequence(row['full_molecule_sequence'], sequencing_type)
            quality_scores = "I" * len(sequence)  # Dummy quality scores

            f.write(f"{read_name}\n")
            f.write(f"{sequence}\n")
            f.write(f"+\n")
            f.write(f"{quality_scores}\n")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="convert dataframe to FASTQ file.")
    parser.add_argument("--input_df", help="full path to the input dataframe")    
    parser.add_argument("--seq_type", choices=["RNA", "cDNA"], default='RNA',help="Type of sequencing: RNA or cDNA")
    parser.add_argument('--seq_depth', type=int, default=20000000, help='total library sequencing depth')

    parser.add_argument('--tpm_lower_limit', type=int, default=5, help='lower possible CPM per gene')
    parser.add_argument('--tpm_upper_limit', type=int, default=200, help='higher possible CPM per gene')

    parser.add_argument('-o', type=str, default='./', help='output path')
    parser.add_argument("--seed", type=int, help="random seed for reproducibility")

    args = parser.parse_args()

    # Set random seed if specified
    if args.seed is not None:
        random.seed(args.seed)


    # Generate output name and file for output file with read and substitution info
    filename = os.path.splitext(os.path.basename(args.input_df))[0]
    # Extract only the gene ID (first part before the first underscore)
    base_filename = filename.split(".")[0]
    # Construct the output filename
    output_filename = str(args.o) + '/reads/' +base_filename +'_'+args.seq_type +".fastq.gz"

    df = pd.read_csv(args.input_df, delimiter="\t")

    path_to_BGmRNAs=os.path.join(os.path.dirname(args.o), "mRNA", base_filename + "_background.tsv.gz")
    if os.path.exists(path_to_BGmRNAs) and os.path.getsize(path_to_BGmRNAs) > 0:
        df_bg = pd.read_csv(path_to_BGmRNAs delimiter="\t")
        df = pd.concat([df, df_bg], ignore_index=True)



    #Get a random CPM for the gene
    cpm_gene = random.randint(args.tpm_lower_limit, args.tpm_upper_limit)

    scaled_cpm = args.seq_depth * cpm_gene / 10**6
    scaled_cpm=round(scaled_cpm)

    # Sample based on scaling
    if len(df) >= scaled_cpm:
        sampled_df = df.sample(n=scaled_cpm, replace=False)
    else:
        extra_needed=scaled_cpm-len(df)
        additional_samples = df.sample(n=extra_needed, replace=True)
        sampled_df = pd.concat([df, additional_samples], ignore_index=True)



    convert_to_fastq(sampled_df, output_filename, args.seq_type,base_filename)

















