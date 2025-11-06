import pandas as pd
import ast
import random
import argparse
import math
import os
import string
import re
import numpy as np


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="convert dataframe to FASTQ file.")
    parser.add_argument("--input_df", help="full path to the mRNA dataframe")
    parser.add_argument("--o", type=str, default="./", help="output path")
    parser.add_argument("--seed", type=int, help="random seed for reproducibility")

    args = parser.parse_args()
    np.random.seed(args.seed) if args.seed else None

    df = pd.read_csv(args.input_df, sep="\t", comment="#")

    # If incorporated_positions is stored as string, convert to list
    df["incorporated_positions"] = df["incorporated_positions"].apply(
        lambda x: ast.literal_eval(x) if isinstance(x, str) else x
    )

    # Get insert size per molecule
    #insert_size_min, insert_size_max = map(float, args.insert_size.split(","))
    #df["insert_size"] = np.random.randint(insert_size_min, insert_size_max + 1, size=len(df))

    # Keep only rows where incorporated_positions is not empty
    df = df[df["incorporated_positions"].apply(lambda x: len(x) > 0)]

    # Now safely compute the min
    df["min_incorporated"] = df["incorporated_positions"].apply(min)

    # Get the sequence up to the first incorporated position
    df["sequence_up_to_min"] = df.apply(
    	lambda row: row["full_molecule_sequence"][: row["min_incorporated"]],
    	axis=1)

    # Generate the output filename
    filename = os.path.splitext(os.path.basename(args.input_df))[0]
    base_filename = filename.split(".")[0]
    output_filename = os.path.join(args.o, base_filename + "_proseq.tsv.gz")

    # Convert to DataFrame with column names
    df_for_export = pd.DataFrame(
        df,
        columns=[
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
            "min_incorporated",
            "sequence_up_to_min",
        ],
    )

    # Save to a tab-separated file
    df_for_export.to_csv(output_filename, sep="\t", index=False, compression="gzip")















