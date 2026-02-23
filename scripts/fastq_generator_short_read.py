import pandas as pd
import gzip
import random
import argparse
import os
import subprocess
import string
import ast
import numpy as np

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))

def run_cmd(cmd):
    # print(f"[Running] {cmd}") # Uncomment for verbose logging
    subprocess.run(cmd, shell=True, check=True)

def make_random_suffix(length=5):
    return ''.join(random.choices(string.ascii_lowercase, k=length))

def reverse_complement(seq):
    complement = str.maketrans('ACGTacgt', 'TGCAtgca')
    return seq.translate(complement)[::-1]

def process_sequence(sequence, reverse_comp=False):
    return reverse_complement(sequence) if reverse_comp else sequence

def get_absolute_coords(rates_df, read_start, read_end):
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

def process_and_write_fastq(
    df_sampled, 
    lookup_dict, 
    output_prefix, 
    sequencing_type, 
    strandedness, 
    read_length, 
    rates_df, 
    ttseq=False
):
    if sequencing_type == "SE":
        f_out = gzip.open(f"{output_prefix}.fastq.gz", 'wt')
    else:
        f1_out = gzip.open(f"{output_prefix}_R1.fastq.gz", 'wt')
        f2_out = gzip.open(f"{output_prefix}_R2.fastq.gz", 'wt')

    try:
        for row in df_sampled.itertuples():
            ref_data = lookup_dict.get(row.transcript_id)
            if not ref_data: continue 
            
            full_seq = ref_data['seq']
            conv_pos = ref_data['conv']
            inc_pos = ref_data['inc']
            mol_id = ref_data['mol_id']

            try:
                r_start, r_end = map(int, row.read_coordinate_split.split('-'))
            except AttributeError:
                continue

            upstream_start = r_start
            upstream_end = min(r_end, r_start + read_length)
            
            downstream_start = max(r_start, r_end - read_length)
            downstream_end = r_end
            
            inc_upstream = [p for p in inc_pos if upstream_start <= p <= upstream_end]
            inc_downstream = [p for p in inc_pos if downstream_start <= p <= downstream_end]
            
            if ttseq and (len(inc_upstream) == 0 and len(inc_downstream) == 0):
                continue

            seq_upstream = full_seq[upstream_start:upstream_end]
            seq_downstream = full_seq[downstream_start:downstream_end]

            conv_upstream = [p - upstream_start for p in conv_pos if upstream_start <= p <= upstream_end]
            conv_downstream = [p - downstream_start for p in conv_pos if downstream_start <= p <= downstream_end]
            
            inc_upstream_rel = [p - upstream_start for p in inc_upstream]
            inc_downstream_rel = [p - downstream_start for p in inc_downstream]

            random_suffix = make_random_suffix()

            if sequencing_type == "SE":
                target_read = None
                if strandedness == "rf":
                    target_read = "downstream"
                    rev = True
                elif strandedness == "fr":
                    target_read = "upstream"
                    rev = False
                elif strandedness == "unstranded":
                    target_read = "downstream" if random.choice([True, False]) else "upstream"
                    rev = target_read == "downstream"
                
                if target_read == "upstream":
                    final_seq = seq_upstream
                    n_inc = len(inc_upstream_rel)
                    p_inc = inc_upstream_rel
                    n_conv = len(conv_upstream)
                    p_conv = conv_upstream
                    abs_coords = get_absolute_coords(rates_df, upstream_start, upstream_end - 1)
                else:
                    final_seq = seq_downstream
                    n_inc = len(inc_downstream_rel)
                    p_inc = inc_downstream_rel
                    n_conv = len(conv_downstream)
                    p_conv = conv_downstream
                    abs_coords = get_absolute_coords(rates_df, downstream_start, downstream_end - 1)

                final_seq_proc = process_sequence(final_seq, reverse_comp=rev)
                
                conv_str = ','.join(str(pos + 1) for pos in p_conv)
                inc_str = ','.join(str(pos + 1) for pos in p_inc)

                read_name = (
                    f"{mol_id}{random_suffix}_{abs_coords}"
                    f"_ninc{n_inc}:{inc_str}_nsubs{n_conv}:{conv_str}"
                )
                q_scores = "I" * len(final_seq_proc)
                f_out.write(f"{read_name}\n{final_seq_proc}\n+\n{q_scores}\n")

            elif sequencing_type == "PE":
                if strandedness == "rf":
                    s1, rev1 = seq_downstream, True
                    s2, rev2 = seq_upstream, False
                    coords1 = get_absolute_coords(rates_df, downstream_start, downstream_end - 1)
                    coords2 = get_absolute_coords(rates_df, upstream_start, upstream_end - 1)
                    meta1 = (len(inc_downstream_rel), inc_downstream_rel, len(conv_downstream), conv_downstream)
                    meta2 = (len(inc_upstream_rel), inc_upstream_rel, len(conv_upstream), conv_upstream)
                elif strandedness == "fr":
                    s1, rev1 = seq_upstream, False
                    s2, rev2 = seq_downstream, True
                    coords1 = get_absolute_coords(rates_df, upstream_start, upstream_end - 1)
                    coords2 = get_absolute_coords(rates_df, downstream_start, downstream_end - 1)
                    meta1 = (len(inc_upstream_rel), inc_upstream_rel, len(conv_upstream), conv_upstream)
                    meta2 = (len(inc_downstream_rel), inc_downstream_rel, len(conv_downstream), conv_downstream)
                elif strandedness == "unstranded":
                    if random.choice([True, False]):
                        s1, rev1 = seq_downstream, True
                        s2, rev2 = seq_upstream, False
                        coords1 = get_absolute_coords(rates_df, downstream_start, downstream_end - 1)
                        coords2 = get_absolute_coords(rates_df, upstream_start, upstream_end - 1)
                        meta1 = (len(inc_downstream_rel), inc_downstream_rel, len(conv_downstream), conv_downstream)
                        meta2 = (len(inc_upstream_rel), inc_upstream_rel, len(conv_upstream), conv_upstream)
                    else:
                        s1, rev1 = seq_upstream, False
                        s2, rev2 = seq_downstream, True
                        coords1 = get_absolute_coords(rates_df, upstream_start, upstream_end - 1)
                        coords2 = get_absolute_coords(rates_df, downstream_start, downstream_end - 1)
                        meta1 = (len(inc_upstream_rel), inc_upstream_rel, len(conv_upstream), conv_upstream)
                        meta2 = (len(inc_downstream_rel), inc_downstream_rel, len(conv_downstream), conv_downstream)

                seq_r1 = process_sequence(s1, reverse_comp=rev1)
                seq_r2 = process_sequence(s2, reverse_comp=rev2)

                name_r1 = f"{mol_id}{random_suffix}_{coords1}_ninc{meta1[0]}:{','.join(str(p+1) for p in meta1[1])}_nsubs{meta1[2]}:{','.join(str(p+1) for p in meta1[3])}"
                name_r2 = f"{mol_id}{random_suffix}_{coords2}_ninc{meta2[0]}:{','.join(str(p+1) for p in meta2[1])}_nsubs{meta2[2]}:{','.join(str(p+1) for p in meta2[3])}"

                f1_out.write(f"{name_r1}\n{seq_r1}\n+\n{'I'*len(seq_r1)}\n")
                f2_out.write(f"{name_r2}\n{seq_r2}\n+\n{'I'*len(seq_r2)}\n")

    finally:
        if sequencing_type == "SE":
            f_out.close()
        else:
            f1_out.close()
            f2_out.close()

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
    parser.add_argument("--experiment_time", type=int, default=15, help="experiment time in minutes")
    parser.add_argument("--bkg_molecules", type=float, default=0.0, help="Proportion of total reads that should be background (0.0 to 1.0). Default 0.")
    parser.add_argument("--no_sizeselection", action="store_true", help="If specified, do not filter fragments by size")
    parser.add_argument("--no_fragmentation", action="store_true", help="If specified, do not fragment; keep full molecule length")

    args = parser.parse_args()

    if args.seed is not None:
        random.seed(args.seed)
        np.random.seed(args.seed)

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
        cmd_chop += ['--seed', str(args.seed)]
    if args.fragments:
        cmd_chop += ['--fragments with_ground_truth']
    if args.no_sizeselection:
        cmd_chop += ['--no_sizeselection']
    if args.no_fragmentation:
        cmd_chop += ['--no_fragmentation']

    run_cmd(" ".join(cmd_chop))
    
    filename = os.path.splitext(os.path.basename(args.input_df))[0]
    base_filename = filename.split(".")[0]
    output_prefix = f"{args.o}/reads/{base_filename}"
    chopped_coordinates_file_path = f"{args.o}/temp/mRNAs_with_fragments/{base_filename}_fragments.tsv"
    rates_for_gene = f"{args.o}/rate_per_gene/{base_filename}_RatesandTraversalTimes.gtf"

    print("[Python] Loading Gene references and fragments...")
    
    df_ref_raw = pd.read_csv(args.input_df, delimiter="\t")
    for col in ['converted_positions', 'incorporated_positions']:
        df_ref_raw[col] = df_ref_raw[col].apply(lambda x: ast.literal_eval(x) if isinstance(x, str) else (x if isinstance(x, list) else []))
    
    df_ref_raw['transcript_id'] = range(1, len(df_ref_raw) + 1)
    
    clean_ref_dict = {}
    for _, row in df_ref_raw.iterrows():
        clean_ref_dict[row['transcript_id']] = {
            'seq': row['full_molecule_sequence'],
            'conv': row['converted_positions'],
            'inc': row['incorporated_positions'],
            'mol_id': row['molecule_id']
        }
    del df_ref_raw 

    try:
        df_frags = pd.read_csv(chopped_coordinates_file_path, delimiter="\t")
        df_frags['read_coordinates'] = df_frags['read_coordinates'].astype(str).str.split(',')
        df_exploded = df_frags.explode('read_coordinates').reset_index(drop=True)
        df_exploded.rename(columns={'read_coordinates': 'read_coordinate_split'}, inplace=True)
    except (FileNotFoundError, pd.errors.EmptyDataError):
        df_exploded = pd.DataFrame(columns=['transcript_id', 'read_coordinate_split'])

    gtf_path = f"{args.o}/gtf/{base_filename}.tsv.gz"
    if os.path.exists(gtf_path):
        df_gtf = pd.read_csv(gtf_path, delimiter="\t")
        exon_df = df_gtf[df_gtf['feature'] == 'exon']
        gene_length = exon_df["sequence"].str.len().sum() / 1000
    else:
        gene_length = 1.0 

    gene_tpm = np.random.uniform(low=int(args.tpm_lower_limit), high=int(args.tpm_upper_limit))
    seq_depth_million = args.seq_depth / 1e6
    
    total_reads_budget = int(gene_length * gene_tpm * seq_depth_million)
    
    if args.bkg_molecules > 0.0:
        reads_to_get_bg = int(total_reads_budget * args.bkg_molecules)
        reads_to_get_gene = total_reads_budget - reads_to_get_bg
        if reads_to_get_gene < 0: 
            reads_to_get_gene = 0
            reads_to_get_bg = total_reads_budget
    else:
        reads_to_get_bg = 0
        reads_to_get_gene = total_reads_budget
    
    print(f"[Python] Total Target: {total_reads_budget} | Gene Reads: {reads_to_get_gene} | BG Reads: {reads_to_get_bg}")

    if not df_exploded.empty and reads_to_get_gene > 0:
        if len(df_exploded) >= reads_to_get_gene:
            df_gene_final = df_exploded.sample(n=reads_to_get_gene, replace=False)
        else:
            remainder = reads_to_get_gene - len(df_exploded)
            if remainder > 0:
                sampled_remainder = df_exploded.sample(n=remainder, replace=True)
                df_gene_final = pd.concat([df_exploded, sampled_remainder], ignore_index=True)
            else:
                df_gene_final = df_exploded
    else:
        df_gene_final = pd.DataFrame()

    print(f"[Python] Final Gene Reads: {len(df_gene_final)}")

    df_bg_final = pd.DataFrame()
    
    if reads_to_get_bg > 0:
        path_to_BGmRNAs = os.path.join(os.path.dirname(args.o), "mRNA", base_filename + "_background.tsv.gz")
        
        if os.path.exists(path_to_BGmRNAs) and os.path.getsize(path_to_BGmRNAs) > 0:
            
            print(f"[Python] Processing Background. Target BG Reads: {reads_to_get_bg}")

            cmd_chop_bg = [
                f"Rscript {os.path.join(SCRIPT_DIR, 'short_read_chopper.R')}",
                f"--tsv {path_to_BGmRNAs}",
                f"--insert_size {args.insert_size}",
                f"--read_length {args.read_length}",
                f"--threads {args.threads}",
                f"--seq_depth {args.seq_depth}",
                f"--tpm_lower_limit {args.tpm_lower_limit}",
                f"--tpm_upper_limit {args.tpm_upper_limit}",
                f"-o {args.o}"
            ]
            if args.seed:
                cmd_chop_bg += ['--seed', str(args.seed)]
            if args.no_sizeselection:
                cmd_chop_bg += ['--no_sizeselection']
            if args.no_fragmentation:
                cmd_chop_bg += ['--no_fragmentation']
            
            run_cmd(" ".join(cmd_chop_bg))
            
            offset_bg = 10_000_000 
            
            df_bg_ref = pd.read_csv(path_to_BGmRNAs, delimiter="\t")
            for col in ['converted_positions', 'incorporated_positions']:
                df_bg_ref[col] = df_bg_ref[col].apply(lambda x: ast.literal_eval(x) if isinstance(x, str) else [])
            
            df_bg_ref['transcript_id'] = range(1, len(df_bg_ref) + 1)
            for _, row in df_bg_ref.iterrows():
                clean_ref_dict[row['transcript_id'] + offset_bg] = {
                    'seq': row['full_molecule_sequence'],
                    'conv': row['converted_positions'],
                    'inc': row['incorporated_positions'],
                    'mol_id': str(row['molecule_id']) + "_BG_"
                }
            del df_bg_ref

            chopped_bg_path = f"{args.o}/temp/mRNAs_with_fragments/{base_filename}_background_fragments.tsv"
            try:
                df_bg_frags = pd.read_csv(chopped_bg_path, delimiter="\t")
                df_bg_frags['transcript_id'] = df_bg_frags['transcript_id'] + offset_bg
                
                df_bg_frags['read_coordinates'] = df_bg_frags['read_coordinates'].astype(str).str.split(',')
                df_bg_exploded = df_bg_frags.explode('read_coordinates').reset_index(drop=True)
                df_bg_exploded.rename(columns={'read_coordinates': 'read_coordinate_split'}, inplace=True)
                
                print(f"[Python] Raw BG Fragments Available: {len(df_bg_exploded)}")

                if not df_bg_exploded.empty:
                    if len(df_bg_exploded) >= reads_to_get_bg:
                        df_bg_final = df_bg_exploded.sample(n=reads_to_get_bg, replace=False)
                    else:
                        remainder = reads_to_get_bg - len(df_bg_exploded)
                        sampled_remainder = df_bg_exploded.sample(n=remainder, replace=True)
                        df_bg_final = pd.concat([df_bg_exploded, sampled_remainder], ignore_index=True)
                else:
                    print("[Warning] Background R script produced 0 fragments.")

            except (FileNotFoundError, pd.errors.EmptyDataError):
                print("[Warning] Background file empty or missing.")

    df_total = pd.concat([df_gene_final, df_bg_final], ignore_index=True)
    df_total = df_total.sample(frac=1).reset_index(drop=True)

    print(f"[Python] Writing {len(df_total)} reads ({len(df_gene_final)} Gene, {len(df_bg_final)} BG) to FASTQ...")
    
    if os.path.exists(rates_for_gene):
        rates_df = pd.read_csv(rates_for_gene, delimiter="\t")
    else:
        rates_df = pd.DataFrame(columns=['nucleotide_coord', 'chromosome', 'absolute_position', 'strand'])

    process_and_write_fastq(
        df_total, 
        clean_ref_dict, 
        output_prefix, 
        args.seq_type, 
        args.s, 
        args.read_length, 
        rates_df, 
        ttseq=False
    )

    try:
        if os.path.exists(chopped_coordinates_file_path):
            os.remove(chopped_coordinates_file_path)
    except OSError: pass
