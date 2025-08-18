import argparse
import subprocess
import random
import os
import glob
import numpy as np
import pandas as pd
from datetime import datetime

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))

def run_cmd(cmd):
    print(f"[Running] {cmd}")
    subprocess.run(cmd, shell=True, check=True)

#def str_to_float_list(s):
#    return [float(x) for x in s.split(',')]

def main():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    group_general_arguments = parser.add_argument_group('general arguments')
    group_tsv_generation = parser.add_argument_group('tsv generation')
    group_rates_per_region = parser.add_argument_group('rates per region')
    group_labeling_strategy = parser.add_argument_group('labeling strategy')
    group_seq_tech = parser.add_argument_group('sequencing strategy')

    # General arguments
    group_general_arguments.add_argument("--mode", required=True, choices=["fullpipeline", "tsvgeneration", "rates_per_region", "labeling_strategy", "seq_tech"], help="pipeline step to run")
    group_general_arguments.add_argument("-o", type=str, default="./", help="output directory")
    group_general_arguments.add_argument("--seed", type=int, help="random seed for reproducibility")


    # Arguments for tsvgeneration
    group_tsv_generation.add_argument("--gtf", type=str, help="genome annotation file")
    group_tsv_generation.add_argument("--genome_fasta", type=str, help="genome fasta file")
    group_tsv_generation.add_argument("--protein_coding_only", action="store_true", help="only use protein-coding genes")
    group_tsv_generation.add_argument("--threads", type=int, default=1, help="Number of threads")
    group_tsv_generation.add_argument("--n_gene_clusters", type=int, default=6, help="number of gene clusters")
    group_tsv_generation.add_argument("--n", type=int, default=120, help="number of genes")

    # Arguments for rates_per_region
    group_rates_per_region.add_argument("--region_size_range", type=str, help="range of region sizes to split the gene into")
    group_rates_per_region.add_argument('--num_regions', type=int, help='desired number of regions to split the gene into')
    group_rates_per_region.add_argument("--elong_rate_range", default="500,5000", help="range of elongation rates")
    group_rates_per_region.add_argument("--pause_occur_probability", default=0, type=float, help="probability of having a pausing event across the isoform")
    group_rates_per_region.add_argument('--pause_time',default="0.1,0.5", type=str, help="length of the pausing event in minutes", required=False)
    group_rates_per_region.add_argument("--flat_rates", action="store_true", help="all nucleotides have the same elongation rate")

    # Arguments for labeling_strategy
    # group_labeling_strategy.add_argument("--labelingmode", type=str, default="monolabel", choices=["monolabel", "kineticbarcoding"], help="Labeling mode to use")
    group_labeling_strategy.add_argument("--labeling_time", type=int, default=15, help="labeling time in minutes")
    group_labeling_strategy.add_argument("--seq_err", default="0.0001,0.0002", help="sequencing error rate")
    group_labeling_strategy.add_argument('--nt_inc_prob',type=str,default='0.09,0.1',help='Comma-separated range of nucleotide incorporation probability')
    group_labeling_strategy.add_argument('--subs_prob',type=str,default='0.95,1',help='Comma-separated range of nucleotide substitution proportions')
    group_labeling_strategy.add_argument("--sub_type", type=str, default="T,T", help="substitution type")
    group_labeling_strategy.add_argument('--bkg_molecules', type=float, default=0,help='proportion of molecules that are derived from non-labeled RNA')
    group_labeling_strategy.add_argument("--drb", action="store_true", help="DRB treatment experiment for transcription synchronization")

    # Arguments for seq_tech
    group_seq_tech.add_argument("--seq_tech", choices=["longread", "shortread"], default="shortread", help="sequencing technology")
    group_seq_tech.add_argument("--seq_type", choices=["RNA", "cDNA", "SE", "PE"], default="SE",help="sequencing type")
    group_seq_tech.add_argument('--seq_depth', type=int, default=20000000, help='total library sequencing depth')
    group_seq_tech.add_argument("--insert_size", type=str, default="200,300", help="length of insert")
    group_seq_tech.add_argument("--tpm", type=str, default="5,200", help="range of TPMs genes can have")
    group_seq_tech.add_argument("--read_length", type=int, default=100, help="length of each read")
    group_seq_tech.add_argument("--s", choices=["rf", "fr", "unstranded"], default="rf", help="library strandedness")

    args = parser.parse_args()
    output_dir = args.o.rstrip("/")

    # Set random seed if specified
    if args.seed is not None:
        np.random.seed(args.seed)
        random.seed(args.seed)

    #both the number and size of regions were specified but they are mutually exclusive
    region_size_range_provided = args.region_size_range is not None
    num_regions_provided = args.num_regions is not None
    
    if region_size_range_provided and num_regions_provided:
        raise ValueError("Specify only one of --region_size_range or --num_regions, not both.")

    # Print errors if non-compatible options are specified
    #Seq type+strand
    if args.seq_tech == "longread":
        if args.seq_type not in ["RNA", "cDNA"]:
            parser.error("When --seq_tech is longread, --seq_type must be either RNA or cDNA.")
            sys.exit(1)
        if args.insert_size != "200,300" or args.read_length != 100 or args.s != "rf":
            print(
            "Warning: With --seq_tech set to 'longread', the options --insert_size, --read_length, and --s are ignored.",file=sys.stderr)
    if args.seq_tech == "shortread":
        if args.seq_type not in ["SE", "PE"]:
            parser.error("When --seq_tech is shortread, --seq_type must be either SE or PE.")
            sys.exit(1)

    #Not a proportion for bkg molecules
    if not (0.0 <= args.bkg_molecules <= 1.0):
        print("Error: --bkg_molecules must be a float between 0 and 1")
        sys.exit(1)

    #Not a probability for pausing event
    if not (0.0 <= args.pause_occur_probability <= 1.0):
        print("Error: --pause_occur_probability must be between 0 and 1")
        sys.exit(1)

    if args.mode in ["fullpipeline", "tsvgeneration"]:
        cmd_seq = [
            f"Rscript {os.path.join(SCRIPT_DIR, 'seq_and_clustering.R')}",
            f"-g {args.gtf}",
            f"-f {args.genome_fasta}",
            f"-o {output_dir}/",
            f"-t {args.threads}",
            f"-c {args.n_gene_clusters}",
            f"-n {args.n}"
        ]
        if args.protein_coding_only:
            cmd_seq.append("--protein_coding_only")
        if args.seed:
            cmd_seq += ['--seed', args.seed]

        np.random.seed(args.seed)
        random.seed(args.seed)
        run_cmd(" ".join(cmd_seq))

    if args.mode in ["fullpipeline", "rates_per_region"]:
        input_dir = f"{output_dir}/gtf"
        tsv_files = sorted(glob.glob(f"{input_dir}/*.tsv.gz"))
        if not tsv_files:
            print(f"No TSV files found in {input_dir}")
            return
        for tsv in tsv_files:
            cmd_rates = [
                f"python {os.path.join(SCRIPT_DIR, 'rates_per_region.py')}",
                f"--tsv {tsv}",
                f"--elong_rate_range {args.elong_rate_range}",
                f"--pause_time {args.pause_time}",
                f"--pause_occur_probability {args.pause_occur_probability}",
                f"--o {output_dir}"
            ]
            if args.flat_rates:
                cmd_rates.append("--flat_rates")
            if args.seed:
                cmd_rates += ['--seed', args.seed]
            #if args.gene_level:
            #    cmd_rates.append("--gene_level")
            if region_size_range_provided:
                cmd_rates += ['--region_size_range', args.region_size_range]
            elif num_regions_provided:
                cmd_rates += ['--num_regions', str(args.num_regions)]
            else:
                raise ValueError("You must provide either --region_size_range or --num_regions")

            run_cmd(" ".join(cmd_rates))

            #Assign random TPMs per gene. Generate random values between 5 and 200
            #raw_values = np.random.randint(5, 201, size=len(region_files))
            #gene_id = [os.path.splitext(os.path.basename(tsv))[0].split('_')[0] for tsv in region_files] #Clean names to make filtering easier later
            #Get the mean elongation rate
            min_elong_rate, max_elong_rate = args.elong_rate_range.split(',')
            min_elong_rate = int(min_elong_rate)
            max_elong_rate = int(max_elong_rate)
            mean_elong_rate= (min_elong_rate+max_elong_rate)/2
            df_mean_elong = pd.DataFrame({'mean_elong_rate': [mean_elong_rate]})
            df_mean_elong.to_csv(f"{output_dir}/rate_per_gene/mean_elong_rate.tsv", sep="\t", index=False)

            #df_weights = pd.DataFrame({'gene_id': gene_id,'tpm': raw_values,'mean_elong_rate': mean_elong_rate})
            #df_weights.to_csv(f"{output_dir}/rate_per_gene/tpm_per_gene.tsv", sep="\t", index=False)


    if args.mode in ["fullpipeline", "labeling_strategy"]:
        dir_region_files=f"{output_dir}/rate_per_gene/"
        region_files = sorted(glob.glob(os.path.join(dir_region_files, "*VariableElongationRateRegions.gtf")))
        # Create mRNA output directory if it doesn't exist
        mRNA_out_dir = os.path.join(output_dir, "mRNA")
        os.makedirs(mRNA_out_dir, exist_ok=True)


        #Store labeling time
        df_mean_elong = pd.read_csv(f"{output_dir}/rate_per_gene/mean_elong_rate.tsv", sep="\t")
        # Add labeling_time column
        df_mean_elong["labeling_time"] = args.labeling_time
        # Write back to the same file
        df_mean_elong.to_csv(f"{output_dir}/rate_per_gene/mean_elong_rate.tsv", sep="\t", index=False)

        for region_file in region_files:
            gene_id = os.path.basename(region_file).split("_")[0]
            if not gene_id.startswith("ENSG"):
                gene_id = "ENSG" + gene_id
            nt_file = os.path.join(dir_region_files, f"{gene_id}_RatesandTraversalTimes.gtf")

            if not os.path.isfile(nt_file):
                print(f"Warning: nt file missing for {gene_id}, skipping.")
                continue

            treatment='no DRB'
            if args.drb:
                treatment='DRB'

            # if args.labelingmode == "monolabel":
            cmd_label = [
                f"python {os.path.join(SCRIPT_DIR, 'mRNA_generator_monolabeled.py')}",
                f"--region_file {region_file}",
                f"--nt_file {nt_file}",
                f"--l {args.labeling_time}",
                f"--seq_err {args.seq_err}",
                f"--nt_inc_rate {args.nt_inc_prob}",
                f"--subs_rate {args.subs_prob}",
                f"--sub_type {args.sub_type}",
                #f"--seq_depth {args.seq_depth}",
                f"--bkg_molecules {args.bkg_molecules}",
                #f"--path_to_tpm {output_dir}/rate_per_gene/tpm_per_gene.tsv",
                f"--o {mRNA_out_dir}/",
                f"--treatment '{treatment}'"
                ]
            if args.seed:
                cmd_label += ['--seed', args.seed]

            run_cmd(" ".join(cmd_label))

    if args.mode in ["fullpipeline", "seq_tech"]:
        tsv_dir = os.path.join(output_dir, "mRNA")
        tsv_files = sorted(glob.glob(os.path.join(tsv_dir, "ENSG*.tsv.gz")))
        reads_out_dir = os.path.join(output_dir, "reads")
        os.makedirs(reads_out_dir, exist_ok=True)
        tpm_values = args.tpm.split(",")
        tpm_lower_lim = int(tpm_values[0])
        tpm_upper_lim = int(tpm_values[1])

        sub_rate_percent_range = list(map(lambda x: float(x), args.tpm.split(',')))
        libraries_out_dir = os.path.join(output_dir, "final_libraries")
        os.makedirs(libraries_out_dir, exist_ok=True)
        today_str = datetime.today().strftime("%Y%m%d")

        for tsv in tsv_files:
            gene_id = os.path.basename(tsv).split("_")[0]
            print(f"Running seq_tech for gene: {gene_id}")

            if args.seq_tech == "longread":
                cmd_seqtech = [
                    f"python {os.path.join(SCRIPT_DIR, 'fastq_generator_long_read.py')}",
                    f"--input_df {tsv}",
                    f"--seq_type {args.seq_type}",
                    f"--seq_depth {args.seq_depth}",
                    f"--tpm_lower_limit {tpm_lower_lim}",
                    f"--tpm_upper_limit {tpm_upper_lim}",
                    f"-o {reads_out_dir}/"
                ]
                if args.seed:
                    cmd_seqtech += ['--seed', args.seed]

                run_cmd(" ".join(cmd_seqtech))
                output_filename = f"{args.seq_tech}_{args.seq_type}_{today_str}.fastq.gz"
                


            else:
                cmd_seqtech = [
                    f"python {os.path.join(SCRIPT_DIR, 'fastq_generator_short_read.py')}",
                    f"--input_df {tsv}",
                    f"--insert_size {args.insert_size}",
                    f"--read_length {args.read_length}",
                    f"--seq_type {args.seq_type}",
                    f"--threads {args.threads}",
                    f"--seq_depth {args.seq_depth}",
                    f"--tpm_lower_limit {tpm_lower_lim}",
                    f"--tpm_upper_limit {tpm_upper_lim}",
                    f"--s {args.s}",
                    f"--o {output_dir}/"
                ]
                if args.seed:
                    cmd_seqtech += ['--seed', args.seed]

                run_cmd(" ".join(cmd_seqtech))
                #Concatenate libraries
                output_filename_R1 = f"{args.seq_tech}_{args.seq_type}_{today_str}_R1.fastq.gz"
                output_filename_R2 = f"{args.seq_tech}_{args.seq_type}_{today_str}_R2.fastq.gz"

        if args.seq_tech == "longread":
            cmd_catlibraries = [f'cat "{reads_out_dir}"/*_{args.seq_type}.fastq.gz > "{libraries_out_dir}/{output_filename}"']
            run_cmd(" ".join(cmd_catlibraries))
        else:
            cmd_catlibraries_R1 = [f'cat "{reads_out_dir}"/*_R1.fastq.gz > "{libraries_out_dir}/{output_filename_R1}"']
            run_cmd(" ".join(cmd_catlibraries_R1))
            cmd_catlibraries_R2 = [f'cat "{reads_out_dir}"/*_R2.fastq.gz > "{libraries_out_dir}/{output_filename_R2}"']
            run_cmd(" ".join(cmd_catlibraries_R2))
        #Remove reads after concatenation
        cmd_clean_indiv_reads = [f'/bin/rm "{reads_out_dir}"/*fastq']
        run_cmd(" ".join(cmd_clean_indiv_reads))


        

if __name__ == "__main__":
    main()



