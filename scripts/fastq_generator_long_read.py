import pandas as pd
import random
import string
import argparse
import gzip
import os
import numpy as np

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

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="convert dataframe to FASTQ file.")
    parser.add_argument("--input_df", help="full path to the input dataframe")    
    parser.add_argument("--seq_type", choices=["RNA", "cDNA"], default='RNA',help="Type of sequencing: RNA or cDNA")
    parser.add_argument('--seq_depth', type=int, default=20000000, help='total library sequencing depth')

    parser.add_argument('--tpm_lower_limit', type=int, default=5, help='lower possible CPM per gene')
    parser.add_argument('--tpm_upper_limit', type=int, default=200, help='higher possible CPM per gene')

    parser.add_argument('-o', type=str, default='./', help='output path')
    parser.add_argument("--seed", type=int, help="random seed for reproducibility")
    
    # NEW ARGUMENT
    parser.add_argument("--bkg_molecules", type=float, default=0.0, 
                        help="Proportion of total reads that should be background (0.0 to 1.0). Default 0.")

    args = parser.parse_args()

    # Set random seed if specified
    if args.seed is not None:
        random.seed(args.seed)
        np.random.seed(args.seed)

    # Generate output name and file for output file with read and substitution info
    filename = os.path.splitext(os.path.basename(args.input_df))[0]
    base_filename = filename.split(".")[0]
    output_filename = os.path.join(args.o, 'reads', f"{base_filename}_{args.seq_type}.fastq.gz")
    
    # Ensure output directory exists
    os.makedirs(os.path.dirname(output_filename), exist_ok=True)

    # ==========================
    # 1. CALCULATE READ COUNTS
    # ==========================
    # Get a random CPM for the gene
    cpm_gene = random.randint(args.tpm_lower_limit, args.tpm_upper_limit)
    
    # Calculate desired Gene Reads
    reads_to_get_gene = int(round(args.seq_depth * cpm_gene / 10**6))
    
    # Calculate desired Background Reads based on proportion
    reads_to_get_bg = 0
    if args.bkg_molecules > 0.0:
        if args.bkg_molecules >= 1.0:
            reads_to_get_bg = 1000 # Fallback default if only background requested
            reads_to_get_gene = 0
        elif reads_to_get_gene > 0:
            ratio = args.bkg_molecules / (1.0 - args.bkg_molecules)
            reads_to_get_bg = int(reads_to_get_gene * ratio)

    print(f"[LongRead] Target: {reads_to_get_gene} gene reads, {reads_to_get_bg} background reads.")

    # ==========================
    # 2. STREAMING OUTPUT
    # ==========================
    with gzip.open(output_filename, 'wt') as f_out:
        
        # --- PROCESS MAIN GENE ---
        if reads_to_get_gene > 0:
            df = pd.read_csv(args.input_df, delimiter="\t")
            
            # Sampling logic
            if not df.empty:
                if len(df) >= reads_to_get_gene:
                    sampled_df = df.sample(n=reads_to_get_gene, replace=False)
                else:
                    extra_needed = reads_to_get_gene - len(df)
                    additional_samples = df.sample(n=extra_needed, replace=True)
                    sampled_df = pd.concat([df, additional_samples], ignore_index=True)
                
                # Write to FASTQ immediately
                for row in sampled_df.itertuples():
                    molecule_id = id_generator()
                    read_name = f"{molecule_id}_{base_filename}"
                    
                    sequence = process_sequence(row.full_molecule_sequence, args.seq_type)
                    quality_scores = "I" * len(sequence)

                    f_out.write(f"{read_name}\n{sequence}\n+\n{quality_scores}\n")
            del df # Free memory

        # --- PROCESS BACKGROUND ---
        if reads_to_get_bg > 0:
            path_to_BGmRNAs = os.path.join(os.path.dirname(args.o), "mRNA", base_filename + "_background.tsv.gz")
            
            if os.path.exists(path_to_BGmRNAs) and os.path.getsize(path_to_BGmRNAs) > 0:
                try:
                    df_bg = pd.read_csv(path_to_BGmRNAs, delimiter="\t")
                    
                    if not df_bg.empty:
                        # Sampling logic for BG
                        if len(df_bg) >= reads_to_get_bg:
                            sampled_bg = df_bg.sample(n=reads_to_get_bg, replace=False)
                        else:
                            extra_needed = reads_to_get_bg - len(df_bg)
                            additional_samples = df_bg.sample(n=extra_needed, replace=True)
                            sampled_bg = pd.concat([df_bg, additional_samples], ignore_index=True)

                        # Write to FASTQ immediately
                        for row in sampled_bg.itertuples():
                            molecule_id = id_generator()
                            # Append BG tag to name? Optional but helpful
                            read_name = f"{molecule_id}_{base_filename}_BG"
                            
                            sequence = process_sequence(row.full_molecule_sequence, args.seq_type)
                            quality_scores = "I" * len(sequence)

                            f_out.write(f"{read_name}\n{sequence}\n+\n{quality_scores}\n")
                    del df_bg
                except Exception as e:
                    print(f"Warning: Failed to process background file: {e}")
