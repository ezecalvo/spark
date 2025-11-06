import argparse
import subprocess
import random
import os
import glob
import sys
import numpy as np
import pandas as pd
from datetime import datetime
from multiprocessing import Pool

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
# Define args as a global variable, to be assigned in main()
args = None 

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# ALL WORKER FUNCTIONS ARE NOW DEFINED AT THE TOP LEVEL OF THE SCRIPT
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

def run_cmd(cmd):
    """Executes a command and checks for errors."""
    print(f"[Running] {cmd}")
    subprocess.run(cmd, shell=True, check=True)

def process_rates_file(tsv):
    """Worker for the rates_per_region step."""
    cmd_rates = [
        f"python {os.path.join(SCRIPT_DIR, 'rates_per_region.py')}",
        f"--tsv {tsv}", f"--elong_rate_range {args.elong_rate_range}",
        f"--pause_time {args.pause_time}", f"--o {args.o.rstrip('/')}",
        f"--promoter_pause_position {args.promoter_pause_position}",
        f"--promoter_pause_duration {args.promoter_pause_duration}"
    ]
    if args.flat_rates: cmd_rates.append("--flat_rates")
    if args.seed: cmd_rates.extend(['--seed', str(args.seed)])
    if args.region_size_range: cmd_rates.extend(['--region_size_range', args.region_size_range])
    elif args.num_regions: cmd_rates.extend(['--num_regions', str(args.num_regions)])
    else: raise ValueError("You must provide either --region_size_range or --num_regions")
    if args.num_pauses is not None and args.num_pauses > 0:
        cmd_rates.append(f"--num_pauses {args.num_pauses}")
    run_cmd(" ".join(cmd_rates))

def process_mrna_file(region_file):
    """Worker for the mRNAgeneration step."""
    output_dir = args.o.rstrip("/")
    dir_region_files = f"{output_dir}/rate_per_gene/"
    mRNA_out_dir = os.path.join(output_dir, "mRNA")
    
    gene_id = os.path.basename(region_file).split("_")[0]
    if not gene_id.startswith("ENSG"): gene_id = "ENSG" + gene_id
    nt_file = os.path.join(dir_region_files, f"{gene_id}_RatesandTraversalTimes.gtf")

    if not os.path.isfile(nt_file):
        print(f"Warning: nt file missing for {gene_id}, skipping.")
        return

    cmd_label = [
        f"python {os.path.join(SCRIPT_DIR, 'mRNA_generator.py')}",
        f"--region_file {region_file}", f"--nt_file {nt_file}",
        f"--experiment_time {args.experiment_time}", f"--mRNAcount 5000",
        f"--seq_err {args.seq_err}", f"--nt_inc_rate {args.nt_inc_prob}",
        f"--subs_rate {args.subs_prob}", f"--labeling_base {args.labeling_base}",
        f"--bkg_molecules {args.bkg_molecules}", f"--initi_rate {args.initi_rate}",
        f"--o {mRNA_out_dir}/"
    ]
    if args.seed: cmd_label.extend(['--seed', str(args.seed)])
    if args.drb: cmd_label.append('--drb')
    if args.sub_base: cmd_label.extend(['--sub_base', args.sub_base])
    run_cmd(" ".join(cmd_label))

def mrna_to_proseq(tsv):
    """Worker for PRO-seq step 1."""
    proseq_mrnas_dirout = os.path.join(args.o.rstrip('/'), "proseqmRNAs")
    cmd = [f"python {os.path.join(SCRIPT_DIR, 'mRNA_to_proseq.py')}", f"--input_df {tsv}", f"--o {proseq_mrnas_dirout}"]
    if args.seed: cmd.extend(['--seed', str(args.seed)])
    run_cmd(" ".join(cmd))

def proseq_to_fastq(tsv):
    """Worker for PRO-seq step 2."""
    tpm_lower, tpm_upper = map(int, args.tpm.split(","))
    cmd = [f"python {os.path.join(SCRIPT_DIR, 'fastq_generator_proseq.py')}",
           f"--input_df {tsv}", f"--insert_size {args.insert_size}", f"--read_length {args.read_length}",
           f"--seq_type {args.seq_type}", f"--threads 1", f"--seq_depth {args.seq_depth}",
           f"--tpm_lower_limit {tpm_lower}", f"--tpm_upper_limit {tpm_upper}", f"--s {args.s}", f"--o {args.o.rstrip('/')}/", f"--bkg_molecules {args.bkg_molecules}"]
    if args.seed: cmd.extend(['--seed', str(args.seed)])
    run_cmd(" ".join(cmd))

def mrna_to_mnetseq(tsv):
    """Worker for mNET-seq step 1."""
    mnetseq_mrnas_dirout = os.path.join(args.o.rstrip('/'), "mnetseqmRNAs")
    cmd = [f"python {os.path.join(SCRIPT_DIR, 'mRNA_to_mnetseq.py')}", f"--input_df {tsv}", f"--o {mnetseq_mrnas_dirout}"]
    if args.seed: cmd.extend(['--seed', str(args.seed)])
    run_cmd(" ".join(cmd))
    
def mnetseq_to_fastq(tsv):
    """Worker for mNET-seq step 2."""
    tpm_lower, tpm_upper = map(int, args.tpm.split(","))
    cmd = [f"python {os.path.join(SCRIPT_DIR, 'fastq_generator_mnetseq.py')}",
           f"--input_df {tsv}", f"--insert_size {args.insert_size}", f"--read_length {args.read_length}",
           f"--seq_type {args.seq_type}", f"--threads 1", f"--seq_depth {args.seq_depth}",
           f"--tpm_lower_limit {tpm_lower}", f"--tpm_upper_limit {tpm_upper}", f"--s {args.s}", f"--o {args.o.rstrip('/')}/"]
    if args.seed: cmd.extend(['--seed', str(args.seed)])
    run_cmd(" ".join(cmd))

def nascent_to_fastq(tsv):
    """Worker for nascentRNA-pd."""
    output_dir = args.o.rstrip('/')
    tpm_lower, tpm_upper = map(int, args.tpm.split(","))
    if args.seq_tech == "longread":
        cmd = [f"python {os.path.join(SCRIPT_DIR, 'fastq_generator_long_read.py')}",
               f"--input_df {tsv}", f"--seq_type {args.seq_type}", f"--seq_depth {args.seq_depth}",
               f"--tpm_lower_limit {tpm_lower}", f"--tpm_upper_limit {tpm_upper}", f"-o {output_dir}"]
    else: # shortread
        cmd = [f"python {os.path.join(SCRIPT_DIR, 'fastq_generator_short_read.py')}",
               f"--input_df {tsv}", f"--insert_size {args.insert_size}", f"--read_length {args.read_length}",
               f"--seq_type {args.seq_type}", f"--threads 1", f"--seq_depth {args.seq_depth}",
               f"--tpm_lower_limit {tpm_lower}", f"--tpm_upper_limit {tpm_upper}", f"--s {args.s}", f"--o {output_dir}/"]
        if args.fragments: cmd.append('--fragments')
    if args.seed: cmd.extend(['--seed', str(args.seed)])

    run_cmd(" ".join(cmd))

def ttseq_to_fastq(tsv):
    """Worker for TT-seq."""
    tpm_lower, tpm_upper = map(int, args.tpm.split(","))
    cmd = [f"python {os.path.join(SCRIPT_DIR, 'mRNA_to_reads_ttseq.py')}",
           f"--input_df {tsv}", f"--insert_size {args.insert_size}", f"--read_length {args.read_length}",
           f"--seq_type {args.seq_type}", f"--threads 1", f"--seq_depth {args.seq_depth}",
           f"--tpm_lower_limit {tpm_lower}", f"--tpm_upper_limit {tpm_upper}", f"--s {args.s}", f"--o {args.o.rstrip('/')}/", f"--bkg_molecules {args.bkg_molecules}"]
    if args.fragments: cmd.append('--fragments')
    if args.seed: cmd.extend(['--seed', str(args.seed)])
    run_cmd(" ".join(cmd))


def main():
    global args # Declare that this function will assign to the global 'args'
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    # (Argument parsing groups remain unchanged)
    group_general_arguments = parser.add_argument_group('general arguments')
    group_tsv_generation = parser.add_argument_group('tsv generation')
    group_rates_per_region = parser.add_argument_group('rates per region')
    group_mRNAgeneration = parser.add_argument_group('mRNA generation')
    group_seq_tech = parser.add_argument_group('sequencing strategy')

    group_general_arguments.add_argument("--mode", required=True, choices=["fullpipeline", "tsvgeneration", "rates_per_region", "mRNAgeneration", "seq_tech"], help="pipeline step to run")
    group_general_arguments.add_argument("-o", type=str, default="./", help="output directory")
    group_general_arguments.add_argument("--seed", type=int, help="random seed for reproducibility")
    group_general_arguments.add_argument("--experiment_type", required=True, choices=["nascentrnapd", "proseq","ttseq","mnetseq"], help="type of experiment to simulate")
    group_general_arguments.add_argument("--threads", type=int, default=1, help="Number of threads for parallel processing")

    group_tsv_generation.add_argument("--gtf", type=str, help="genome annotation file")
    group_tsv_generation.add_argument("--genome_fasta", type=str, help="genome fasta file")
    group_tsv_generation.add_argument("--protein_coding_only", action="store_true", help="only use protein-coding genes")
    group_tsv_generation.add_argument("--n_gene_clusters", type=int, default=6, help="number of gene clusters")
    group_tsv_generation.add_argument("--n", type=int, default=120, help="number of genes")

    group_rates_per_region.add_argument("--region_size_range", type=str, help="range of region sizes to split the gene into")
    group_rates_per_region.add_argument('--num_regions', type=int, help='desired number of regions to split the gene into')
    group_rates_per_region.add_argument("--elong_rate_range", default="500,5000", help="range of elongation rates")
    group_rates_per_region.add_argument('--num_pauses', type=int, help='number of pauses across the gene',default=0)
    group_rates_per_region.add_argument('--pause_time',default="0.1,0.5", type=str, help="length of the pausing event in minutes", required=False)
    group_rates_per_region.add_argument("--flat_rates", action="store_true", help="all nucleotides have the same elongation rate")
    group_rates_per_region.add_argument('--promoter_pause_position', type=str, default='30,50', help='range of positions for TSS pausing')
    group_rates_per_region.add_argument('--promoter_pause_duration', type=str, default='0,0', help='length of promoter pausing in minutes')
    group_rates_per_region.add_argument('--initi_rate',type=str,default='10,60',help='comma-separated range of seconds for initiation rate')

    group_mRNAgeneration.add_argument("--experiment_time", type=int, default=15, help="experiment time in minutes")
    group_mRNAgeneration.add_argument("--seq_err", default="0.0001,0.0002", help="sequencing error rate")
    group_mRNAgeneration.add_argument('--nt_inc_prob',type=str,default='0.09,0.1',help='comma-separated range of nucleotide incorporation probability')
    group_mRNAgeneration.add_argument('--subs_prob',type=str,default='0.95,1',help='comma-separated range of nucleotide substitution probability if simulating nucleotide recoding')
    group_mRNAgeneration.add_argument('--labeling_base', type=str, default='T', choices=["A", "C", "T", "G"],help='Specify the base analog that will get incorporated into the nascentRNA', required=False)
    group_mRNAgeneration.add_argument('--sub_base', type=str, choices=["A", "C", "T", "G"], help='Specify the identity of the base after conversion', required=False)
    group_mRNAgeneration.add_argument('--bkg_molecules', type=float, default=0,help='proportion of mRNA molecules that are derived from non-labeled RNA')
    group_mRNAgeneration.add_argument("--drb", action="store_true", help="DRB treatment experiment for transcription synchronization")

    group_seq_tech.add_argument("--seq_tech", choices=["longread", "shortread"], default="shortread", help="sequencing technology")
    group_seq_tech.add_argument("--seq_type", choices=["RNA", "cDNA", "SE", "PE"], default="SE",help="sequencing type")
    group_seq_tech.add_argument('--seq_depth', type=int, default=20000000, help='total library sequencing depth')
    group_seq_tech.add_argument("--insert_size", type=str, default="200,300", help="length of insert")
    group_seq_tech.add_argument("--tpm", type=str, default="5,200", help="range of TPMs genes can have")
    group_seq_tech.add_argument("--read_length", type=int, default=100, help="length of each read")
    group_seq_tech.add_argument("--s", choices=["rf", "fr", "unstranded"], default="rf", help="library strandedness")
    group_seq_tech.add_argument("--fragments", action="store_true", help="export a ground truth for fragmentation & size selection")

    args = parser.parse_args()
    output_dir = args.o.rstrip("/")

    if args.seed is not None:
        np.random.seed(args.seed)
        random.seed(args.seed)
    
    # (Input validation checks remain unchanged)
    if args.region_size_range and args.num_regions:
        raise ValueError("Specify only one of --region_size_range or --num_regions, not both.")
    if args.n_gene_clusters > args.n:
        sys.exit(f"Error: --n_gene_clusters ({args.n_gene_clusters}) cannot be greater than -n ({args.n})")
    if args.experiment_type == 'mnetseq':
        print("Warning: nt_inc_prob, subs_prob, or sub_base are not used for mNET-seq simulations")
    if args.seq_tech == "longread" and args.seq_type not in ["RNA", "cDNA"]:
        parser.error("When --seq_tech is longread, --seq_type must be either RNA or cDNA.")
    if args.seq_tech == "shortread" and args.seq_type not in ["SE", "PE"]:
        parser.error("When --seq_tech is shortread, --seq_type must be either SE or PE.")
    if not (0.0 <= args.bkg_molecules <= 1.0):
        sys.exit("Error: --bkg_molecules must be a float between 0 and 1")


    # --- Pipeline Steps ---
    if args.mode in ["fullpipeline", "tsvgeneration"]:
        cmd_seq = [
            f"Rscript {os.path.join(SCRIPT_DIR, 'seq_and_clustering.R')}",
            f"-g {args.gtf}", f"-f {args.genome_fasta}", f"-o {output_dir}/",
            f"-t {args.threads}", f"-c {args.n_gene_clusters}", f"-n {args.n}"
        ]
        if args.protein_coding_only: cmd_seq.append("--protein_coding_only")
        if args.seed: cmd_seq.extend(['--seed', str(args.seed)])
        run_cmd(" ".join(cmd_seq))

    if args.mode in ["fullpipeline", "rates_per_region"]:
        input_dir = f"{output_dir}/gtf"
        tsv_files = sorted(glob.glob(f"{input_dir}/*.tsv.gz"))
        if not tsv_files: sys.exit(f"No TSV files found in {input_dir}")
        with Pool(processes=args.threads) as pool:
            pool.map(process_rates_file, tsv_files)

    if args.mode in ["fullpipeline", "mRNAgeneration"]:
        dir_region_files=f"{output_dir}/rate_per_gene/"
        region_files = sorted(glob.glob(os.path.join(dir_region_files, "*VariableElongationRateRegions.gtf")))
        mRNA_out_dir = os.path.join(output_dir, "mRNA")
        os.makedirs(mRNA_out_dir, exist_ok=True)
        with Pool(processes=args.threads) as pool:
            pool.map(process_mrna_file, region_files)
    
    if args.mode in ["fullpipeline", "seq_tech"]:
        today_str = datetime.today().strftime("%Y%m%d")
        reads_out_dir = os.path.join(output_dir, "reads")
        libraries_out_dir = os.path.join(output_dir, "final_libraries")
        os.makedirs(reads_out_dir, exist_ok=True)
        os.makedirs(libraries_out_dir, exist_ok=True)

    # --- Sequencing Technology Simulation (Parallelized) ---
    if args.experiment_type == 'proseq' and args.mode in ["fullpipeline", "seq_tech"]:
        proseq_mrnas_dirout = os.path.join(output_dir, "proseqmRNAs")
        os.makedirs(proseq_mrnas_dirout, exist_ok=True)
        tsv_files = [f for f in sorted(glob.glob(os.path.join(output_dir, "mRNA", "ENSG*.tsv.gz"))) if "background" not in os.path.basename(f)]
        with Pool(args.threads) as pool:
            pool.map(mrna_to_proseq, tsv_files)
        
        proseq_files = [f for f in sorted(glob.glob(os.path.join(proseq_mrnas_dirout, "ENSG*_proseq.tsv.gz"))) if "background" not in os.path.basename(f)]
        with Pool(args.threads) as pool:
            pool.map(proseq_to_fastq, proseq_files)

    elif args.experiment_type == 'mnetseq' and args.mode in ["fullpipeline", "seq_tech"]:
        mnetseq_mrnas_dirout = os.path.join(output_dir, "mnetseqmRNAs")
        os.makedirs(mnetseq_mrnas_dirout, exist_ok=True)
        tsv_files = [f for f in sorted(glob.glob(os.path.join(output_dir, "mRNA", "ENSG*.tsv.gz"))) if "background" not in os.path.basename(f)]
        with Pool(args.threads) as pool:
            pool.map(mrna_to_mnetseq, tsv_files)

        mnetseq_files = [f for f in sorted(glob.glob(os.path.join(mnetseq_mrnas_dirout, "ENSG*_mnetseq.tsv.gz"))) if "background" not in os.path.basename(f)]
        with Pool(args.threads) as pool:
            pool.map(mnetseq_to_fastq, mnetseq_files)

    elif args.experiment_type == 'nascentrnapd' and args.mode in ["fullpipeline", "seq_tech"]:
        tsv_files = [f for f in sorted(glob.glob(os.path.join(output_dir, "mRNA", "ENSG*.tsv.gz"))) if "background" not in os.path.basename(f)]
        with Pool(args.threads) as pool:
            pool.map(nascent_to_fastq, tsv_files)

    elif args.experiment_type == 'ttseq' and args.mode in ["fullpipeline", "seq_tech"]:
        tsv_files = [f for f in sorted(glob.glob(os.path.join(output_dir, "mRNA", "ENSG*.tsv.gz"))) if "background" not in os.path.basename(f)]
        with Pool(args.threads) as pool:
            pool.map(ttseq_to_fastq, tsv_files)
            
    # --- Final Library Concatenation ---
    if args.mode in ["fullpipeline", "seq_tech"]:
        if not glob.glob(f"{reads_out_dir}/*fastq*"):
            print("No read files found to concatenate. Skipping final step.")
            return

        if args.seq_tech == "longread":
            output_filename = f"{args.seq_tech}_{args.experiment_type}_{args.seq_type}_{today_str}.fastq.gz"
            run_cmd(f'cat "{reads_out_dir}"/*.fastq.gz > "{libraries_out_dir}/{output_filename}"')
        elif args.seq_tech == "shortread" and args.seq_type == "PE":
            output_r1 = f"{args.seq_tech}_{args.experiment_type}_{args.seq_type}_{today_str}_R1.fastq.gz"
            output_r2 = f"{args.seq_tech}_{args.experiment_type}_{args.seq_type}_{today_str}_R2.fastq.gz"
            run_cmd(f'cat "{reads_out_dir}"/*_R1.fastq.gz > "{libraries_out_dir}/{output_r1}"')
            run_cmd(f'cat "{reads_out_dir}"/*_R2.fastq.gz > "{libraries_out_dir}/{output_r2}"')
        else: # shortread SE
            output_se = f"{args.seq_tech}_{args.experiment_type}_{args.seq_type}_{today_str}.fastq.gz"
            run_cmd(f'cat "{reads_out_dir}"/*.fastq.gz > "{libraries_out_dir}/{output_se}"')
            
        run_cmd(f'rm -f "{reads_out_dir}"/*fastq*')
        temp_dir = f"{output_dir}/temp"
        if os.path.isdir(temp_dir):
            run_cmd(f'rm -rf "{temp_dir}"')

if __name__ == "__main__":
    main()

