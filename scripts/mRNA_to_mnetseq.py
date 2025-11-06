import pandas as pd
import ast
import random
import argparse
import math
import os
import string
import re
import numpy as np


#Add function to chop the bg molecules if needed. This function will chop in 35-100nt fragments the full length bg molecule
def chop_sequence(seq):
    fragment_length = random.randint(35, 100)
    max_start_pos = len(seq) - fragment_length
    start_pos = random.randint(0, max_start_pos)
    fragment = seq[start_pos:start_pos + fragment_length]
    return [(start_pos, start_pos + fragment_length, fragment)]


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="convert mRNAs to mNET-seq fragments")
    parser.add_argument("--input_df", help="full path to the mRNA dataframe")
    parser.add_argument("--o", type=str, default="./", help="output path")
    parser.add_argument("--seed", type=int, help="random seed for reproducibility")


    args = parser.parse_args()

    # Generate the output filename
    filename = os.path.splitext(os.path.basename(args.input_df))[0]
    base_filename = filename.split(".")[0]
    output_filename = os.path.join(args.o, base_filename + "_mnetseq.tsv.gz")



    
    np.random.seed(args.seed) if args.seed else None

    df = pd.read_csv(args.input_df, sep="\t", comment="#")


    #Get max length for the transcript to remove mRNAs that are done elongating at this time
    path_to_rates = os.path.join(os.path.dirname(args.o), "rate_per_gene", base_filename + "_RatesandTraversalTimes.gtf")

    rates_df = pd.read_csv(path_to_rates, sep="\t", comment="#")
    rates_df["nucleotide_coord"] = rates_df["nucleotide_coord"].astype(int)
    max_length=rates_df["nucleotide_coord"].max()
    
    df = df[df['stop_label_pos'] < (max_length - 1)]

    #Get the mnase digestion lengths
    random_offsets = np.random.randint(35, 101, size=len(df))

    #Get the position for the fragments
    df["fragment_start_pos"] = (df["stop_label_pos"] - random_offsets).clip(lower=0)
    df["fragment_end_pos"] = df["stop_label_pos"]

    # Extract substring
    df["parsed_sequence"] = [
    seq[start:end] for seq, start, end in zip(
        df["full_molecule_sequence"], df["fragment_start_pos"], df["fragment_end_pos"]
        )
    ]

    #Only keep fragments longer than 10nt
    df = df[df["parsed_sequence"].str.len() > 10]

    # Convert to DataFrame with column names
    df_for_export = pd.DataFrame(
        df,
        columns=[
            "initiation_time",
            "molecule_id",
            "strand",
            "sub_rate_percent_range",
            "fragment_start_pos",
            "fragment_end_pos",
            "parsed_sequence",
        ],
    )


    #Add background molecules if simulated
    path_to_BGmRNAs=os.path.join(os.path.dirname(args.o), "mRNA", base_filename + "_background.tsv.gz")

    if os.path.exists(path_to_BGmRNAs) and os.path.getsize(path_to_BGmRNAs) > 0:
        df_bg = pd.read_csv(path_to_BGmRNAs, sep="\t", comment="#")
        df_bg['molecule_id'] = df_bg['molecule_id'].astype(str) + '_BG_'#Add to the molecule ID that this is a background molecule
        fragments = []
        rows = []
        for idx, row in df_bg.iterrows():
            fragments = chop_sequence(row['full_molecule_sequence'])
            for fstart, fend, fseq in fragments:
                new_row = {
                'initiation_time': row['initiation_time'],
                'molecule_id': row['molecule_id'],
                'strand': row['strand'],
                'sub_rate_percent_range': row['sub_rate_percent_range'],
                'fragment_start_pos': fstart,
                'fragment_end_pos': fend,
                'parsed_sequence': fseq}
            rows.append(new_row)

        parsed_df = pd.DataFrame(rows)
        df_for_export = pd.concat([df_for_export, parsed_df], ignore_index=True)

    # Save to a tab-separated file
    df_for_export.to_csv(output_filename, sep="\t", index=False, compression="gzip")
























