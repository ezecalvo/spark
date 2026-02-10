import pandas as pd
import random
import argparse
import os
import numpy as np
import gzip

def chop_sequence_bg(seq):
    """
    Chops a random fragment (35-100nt) from a background sequence.
    Returns: (start, end, fragment_string)
    """
    seq_len = len(seq)
    fragment_length = random.randint(35, 100)
    
    # If sequence is shorter than fragment length, take the whole thing
    if seq_len <= fragment_length:
        return 0, seq_len, seq

    max_start_pos = seq_len - fragment_length
    start_pos = random.randint(0, max_start_pos)
    end_pos = start_pos + fragment_length
    return start_pos, end_pos, seq[start_pos:end_pos]

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="convert mRNAs to mNET-seq fragments")
    parser.add_argument("--input_df", help="full path to the mRNA dataframe")
    parser.add_argument("--o", type=str, default="./", help="output path")
    parser.add_argument("--seed", type=int, help="random seed for reproducibility")
    # New argument to control background proportion
    parser.add_argument("--bkg_molecules", type=float, default=0.0, 
                        help="Proportion of total reads that should be background (0.0 to 1.0). Default 0.")

    args = parser.parse_args()

    # Generate the output filename
    filename = os.path.splitext(os.path.basename(args.input_df))[0]
    base_filename = filename.split(".")[0]
    output_filename = os.path.join(args.o, base_filename + "_mnetseq.tsv.gz")

    if args.seed:
        np.random.seed(args.seed)
        random.seed(args.seed)

    # 1. Load Rates to get Max Length
    path_to_rates = os.path.join(os.path.dirname(args.o), "rate_per_gene", base_filename + "_RatesandTraversalTimes.gtf")
    
    # Check if rates file exists, otherwise set a safe default or error
    if os.path.exists(path_to_rates):
        rates_df = pd.read_csv(path_to_rates, sep="\t", comment="#")
        rates_df["nucleotide_coord"] = rates_df["nucleotide_coord"].astype(int)
        max_length = rates_df["nucleotide_coord"].max()
    else:
        # Fallback if file missing (though logic implies it should exist)
        max_length = float('inf')

    # Prepare Output File (Write Header)
    output_columns = [
        "initiation_time", "molecule_id", "strand", "sub_rate_percent_range",
        "fragment_start_pos", "fragment_end_pos", "parsed_sequence"
    ]

    with gzip.open(output_filename, 'wt') as f_out:
        f_out.write("\t".join(output_columns) + "\n")

        # ==========================================
        # PART 1: Process Main Gene Molecules
        # ==========================================
        print(f"[mNET-seq] Processing main gene: {base_filename}")
        
        # Load main dataframe
        df = pd.read_csv(args.input_df, sep="\t", comment="#")
        
        # Filter: remove mRNAs that are done elongating
        df = df[df['stop_label_pos'] < (max_length - 1)].copy()

        # Vectorized offset calculation
        # Mnase digestion lengths (35-100)
        random_offsets = np.random.randint(35, 101, size=len(df))

        # Calculate positions
        df["fragment_start_pos"] = (df["stop_label_pos"] - random_offsets).clip(lower=0)
        df["fragment_end_pos"] = df["stop_label_pos"]
        
        # Count valid gene fragments
        gene_fragment_count = 0

        # Iterate and write immediately (Streaming)
        # This prevents storing the generated sequences in memory
        for row in df.itertuples():
            # Extract substring
            # Note: Pandas slicing is inclusive-exclusive, string slicing is inclusive-exclusive
            seq_slice = row.full_molecule_sequence[row.fragment_start_pos : row.fragment_end_pos]
            
            # Only keep fragments longer than 10nt
            if len(seq_slice) > 10:
                gene_fragment_count += 1
                f_out.write(
                    f"{row.initiation_time}\t"
                    f"{row.molecule_id}\t"
                    f"{row.strand}\t"
                    f"{row.sub_rate_percent_range}\t"
                    f"{row.fragment_start_pos}\t"
                    f"{row.fragment_end_pos}\t"
                    f"{seq_slice}\n"
                )
        
        # Free memory of main df
        del df

        # ==========================================
        # PART 2: Process Background Molecules
        # ==========================================
        if args.bkg_molecules > 0.0:
            path_to_BGmRNAs = os.path.join(os.path.dirname(args.o), "mRNA", base_filename + "_background.tsv.gz")
            
            if os.path.exists(path_to_BGmRNAs) and os.path.getsize(path_to_BGmRNAs) > 0:
                
                # Calculate required background count
                # n_bg = n_gene * (P / (1-P))
                if gene_fragment_count > 0 and args.bkg_molecules < 1.0:
                    ratio = args.bkg_molecules / (1.0 - args.bkg_molecules)
                    reads_to_get_bg = int(gene_fragment_count * ratio)
                elif args.bkg_molecules >= 1.0:
                     # If only background requested or no gene reads found
                     reads_to_get_bg = 1000 # Default fallback
                else:
                    reads_to_get_bg = 0

                if reads_to_get_bg > 0:
                    print(f"[mNET-seq] Adding {reads_to_get_bg} background molecules...")
                    
                    df_bg = pd.read_csv(path_to_BGmRNAs, sep="\t", comment="#")
                    
                    # If we need more reads than we have BG molecules, we sample with replacement
                    if len(df_bg) >= reads_to_get_bg:
                        bg_sample = df_bg.sample(n=reads_to_get_bg, replace=False)
                    else:
                        bg_sample = df_bg.sample(n=reads_to_get_bg, replace=True)
                    
                    # Iterate and process
                    for row in bg_sample.itertuples():
                        mol_id = str(row.molecule_id) + '_BG_'
                        
                        fstart, fend, fseq = chop_sequence_bg(row.full_molecule_sequence)
                        
                        f_out.write(
                            f"{row.initiation_time}\t"
                            f"{mol_id}\t"
                            f"{row.strand}\t"
                            f"{row.sub_rate_percent_range}\t"
                            f"{fstart}\t"
                            f"{fend}\t"
                            f"{fseq}\n"
                        )
                    del df_bg
