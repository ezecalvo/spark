import pandas as pd
import ast
import random
import argparse
import os
import numpy as np
import hashlib

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="convert dataframe to PRO-seq compatible coordinates")
    parser.add_argument("--input_df", help="full path to the mRNA dataframe")
    parser.add_argument("--o", type=str, default="./", help="output path")
    parser.add_argument("--seed", type=int, help="random seed for reproducibility")

    args = parser.parse_args()
    
    
    filename = os.path.basename(args.input_df)
    base_filename = filename.replace(".tsv.gz", "")

    if args.seed is not None:
        hasher = hashlib.sha256(base_filename.encode('utf-8'))
        gene_hash = int(hasher.hexdigest(), 16)
        unique_seed = (args.seed + gene_hash) % 4294967295 
        random.seed(unique_seed)
        np.random.seed(unique_seed)

    df = pd.read_csv(args.input_df, sep="\t", comment="#")

    #if incorporated_positions is stored as string, convert to list
    df["incorporated_positions"] = df["incorporated_positions"].apply(
        lambda x: ast.literal_eval(x) if isinstance(x, str) else (x if isinstance(x, list) else [])
    )

    #keep only rows where an analog was actually incorporated. This was my way to code that the polymerase is stalled there
    df = df[df["incorporated_positions"].apply(lambda x: len(x) > 0)].copy()

    #the first position were the nucleotide analog (biotin-NTP) was incorporated is actually the last nucleotide of the RNA
    df["sequence_length"] = df["incorporated_positions"].apply(min)#update seq length for the chopper that comes after tis

    #export
    filename = os.path.splitext(os.path.basename(args.input_df))[0]
    base_filename = filename.split(".")[0]
    output_filename = os.path.join(args.o, base_filename + "_proseq.tsv.gz")

    
    columns_to_keep = [
        "initiation_time",
        "molecule_id",
        "strand",
        "sub_rate_percent_range",
        "start_label_pos",
        "stop_label_pos",
        "converted_positions",
        "incorporated_positions",
        "percentage_seq_err",
        "seq_err_positions",
        "sequence_length"
    ]
    
    final_cols = [c for c in columns_to_keep if c in df.columns]

    df_for_export = df[final_cols]
    df_for_export.to_csv(output_filename, sep="\t", index=False, compression="gzip")
