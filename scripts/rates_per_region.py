import pandas as pd
import random
import argparse
import math
import os

def saveGTF(input_df, output_filename):
    """Saves a DataFrame to a GTF-like tab-separated file."""
    input_df.to_csv(output_filename, sep='\t', index=False)

def getDNAseq(dna_sequence):
    """Concatenates a sequence from a pandas Series."""
    return ''.join(str(i) for i in dna_sequence)

def defineElongRateRegions(input_df_basename, df, dna_sequence,
                           min_length, max_length,
                           min_rate, max_rate,
                           pause_rate,
                           promoter_pause_upstream_limit, promoter_pause_downstream_limit,
                           TSS_pause_duration, num_pauses,
                           num_regions=None):
    """
    Defines regions with varying or flat elongation rates, incorporating random and TSS-specific pauses.
    """
    region_data = []
    gene_length = len(dna_sequence)
    start_coord = 1

    min_rate = float(min_rate)
    max_rate = float(max_rate)

    chromosome = df.iloc[0, 0]
    strand = df.iloc[0, 6]
    abs_start = df.iloc[0, 1]
    abs_end = df.iloc[0, 2]

    # --- Generate region lengths ---
    if num_regions is not None:
        lengths = []
        total = 0
        for i in range(num_regions - 1):
            remaining = gene_length - total - (num_regions - len(lengths) - 1)
            l = random.randint(1, remaining)
            lengths.append(l)
            total += l
        lengths.append(gene_length - total)
    else:
        lengths = []
        total = 0
        while total < gene_length:
            l = random.randint(min_length, max_length)
            if total + l > gene_length:
                l = gene_length - total
            lengths.append(l)
            total += l

    # --- Random pauses ---
    pause_positions = set()
    if num_pauses > 0:
        if num_pauses > gene_length:
            raise ValueError("num_pauses cannot exceed gene length (nt)")
        pause_positions = set(random.sample(range(1, gene_length + 1), num_pauses))

    # --- TSS pause ---
    TSS_pause_position = random.randint(
        promoter_pause_upstream_limit,
        min(promoter_pause_downstream_limit, gene_length)
    )
    all_pauses = set(pause_positions)

    if TSS_pause_duration > 0:
        all_pauses.add(TSS_pause_position)
        TSS_pause_elong_rate = 1 / TSS_pause_duration

    # --- Build regions ---
    rate_initial = round(random.uniform(min_rate, max_rate), 4)
    for region_length in lengths:
        end_coord = start_coord + region_length - 1

        # Pauses inside this region
        pauses_in_region = sorted([p for p in all_pauses if start_coord <= p <= end_coord])

        pos = start_coord
        for pause_pos in pauses_in_region:
            # Normal segment before the pause
            if pos < pause_pos:
                sub_end = pause_pos - 1
                sub_len = sub_end - pos + 1
                sequence = dna_sequence[pos - 1:sub_end]

                # For any pause, the segments before and after have the same flat rate.
                rate_final = rate_initial
                rate_change_per_nt = 0
                time_to_traverse = sub_len / rate_initial if sub_len > 0 and rate_initial > 0 else 0

                region_data.append((chromosome, abs_start, abs_end,
                                    pos, sub_end, strand,
                                    rate_initial, rate_final,
                                    rate_change_per_nt, time_to_traverse,
                                    sequence))
                
                # The rate for the segment after the pause will be the same as the rate before.
                pos = sub_end + 1

            # Pause itself (1 nt)
            sequence = dna_sequence[pause_pos - 1:pause_pos]
            if pause_pos == TSS_pause_position:
                time_to_traverse = 1 / TSS_pause_elong_rate
                rate_val = TSS_pause_elong_rate
            else:
                time_to_traverse = 1 / pause_rate
                rate_val = pause_rate

            region_data.append((chromosome, abs_start, abs_end,
                                pause_pos, pause_pos, strand,
                                rate_val, rate_val, 0,
                                time_to_traverse, sequence))
            pos = pause_pos + 1

        # Final segment after last pause (if any)
        if pos <= end_coord:
            sequence = dna_sequence[pos - 1:end_coord]
            sub_len = end_coord - pos + 1

            # If the whole gene is a single region without pauses, force a flat rate.
            if num_regions == 1 and not all_pauses:
                rate_final = rate_initial
            else:
                rate_final = rate_initial if args.flat_rates else round(random.uniform(min_rate, max_rate), 4)
            
            # Corrected calculation for rate change per nucleotide
            rate_change_per_nt = (rate_final - rate_initial) / (sub_len - 1) if sub_len > 1 else 0
            time_to_traverse = sub_len / ((rate_final + rate_initial) / 2) if sub_len > 0 else 0

            region_data.append((chromosome, abs_start, abs_end,
                                pos, end_coord, strand,
                                rate_initial, rate_final,
                                rate_change_per_nt, time_to_traverse,
                                sequence))
            
            rate_initial = rate_final if not args.flat_rates else round(random.uniform(min_rate, max_rate), 4)

        start_coord = end_coord + 1

    columns = ["chromosome", 'absolute_start', 'absolute_end',
               "region_start_coord", "region_end_coord",
               'strand', "rate_initial", "rate_final",
               "rate_change_per_nt", "time_to_traverse", "sequence"]

    regions_df = pd.DataFrame(region_data, columns=columns)

    # Adjust absolute coordinates based on strand, overwriting original columns to maintain output format
    if not regions_df.empty:
        # Store original gene boundaries before overwriting the columns
        gene_start_col = regions_df['absolute_start'].iloc[0]
        gene_end_col = regions_df['absolute_end'].iloc[0]

        if strand == '+':
            regions_df["absolute_start"] = gene_start_col - 1 + regions_df["region_start_coord"]
            regions_df["absolute_end"]   = gene_start_col - 1 + regions_df["region_end_coord"]
        else: # strand == '-'
            regions_df["absolute_start"] = gene_end_col - regions_df["region_end_coord"] + 1
            regions_df["absolute_end"]   = gene_end_col - regions_df["region_start_coord"] + 1

    return regions_df


def ntTraversalTime(region_specific_eRates_df):
    """Calculates the traversal time for each nucleotide based on region rates."""
    if region_specific_eRates_df.empty:
        return pd.DataFrame()
        
    chromosome = region_specific_eRates_df.iloc[0, 0]
    strand = region_specific_eRates_df.iloc[0, 5]
    abs_pos_start_gene = region_specific_eRates_df.iloc[0, 1] if strand == '+' else region_specific_eRates_df.iloc[0, 2]

    traversal_time_per_nt = []
    cumulative_nt_position = 1

    for region in range(len(region_specific_eRates_df)):
        region_length = (region_specific_eRates_df['region_end_coord'].iloc[region] -
                         region_specific_eRates_df['region_start_coord'].iloc[region] + 1)
        rate_initial = float(region_specific_eRates_df['rate_initial'].iloc[region])
        rate_change_per_nt = float(region_specific_eRates_df['rate_change_per_nt'].iloc[region])

        for nt_position in range(region_length):
            rate_for_this_nt = rate_initial + (rate_change_per_nt * nt_position)
            time_for_this_nt = 1 / rate_for_this_nt if rate_for_this_nt > 0 else float('inf')


            traversal_time_per_nt.append((chromosome, abs_pos_start_gene, strand, region + 1,
                                          cumulative_nt_position, time_for_this_nt,
                                          rate_for_this_nt, rate_change_per_nt))
            cumulative_nt_position += 1

    columns_nt_table = ["chromosome", "absolute_position", 'strand', "region_number", "nucleotide_coord",
                        "time_for_this_nt", "rate_for_this_nt", "rate_change_per_nt"]
    nt_traversal_time_df = pd.DataFrame(traversal_time_per_nt, columns=columns_nt_table)
    
    if not nt_traversal_time_df.empty:
        if strand == '+':
            nt_traversal_time_df["absolute_position"] = (nt_traversal_time_df["absolute_position"].iloc[0] - 1 +
                                                         nt_traversal_time_df["nucleotide_coord"])
        else:
            nt_traversal_time_df["absolute_position"] = (nt_traversal_time_df["absolute_position"].iloc[0] + 1 -
                                                         nt_traversal_time_df["nucleotide_coord"])

    return nt_traversal_time_df

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate variable elongation rate profiles for genes.")
    parser.add_argument('--region_size_range', type=str, help='min,max range of region sizes (in nt)')
    parser.add_argument('--num_regions', type=int, help='desired number of regions to split the gene into')
    parser.add_argument('--elong_rate_range', type=str, default='500,5000', help='range of elongation rates (nt/min)')
    parser.add_argument('--num_pauses', default=0, type=int, help='number of random pause regions across the gene')
    parser.add_argument('--pause_time', type=str, default="0.1,0.5", help='range for length of pausing in minutes')
    parser.add_argument('--tsv', type=str, required=True, help='input TSV file with gene info and sequence')
    parser.add_argument('--o', type=str, default='./', help='output path')
    parser.add_argument('--flat_rates', action='store_true', help='make rates constant within non-pause regions')
    parser.add_argument('--promoter_pause_position', type=str, default='30,50', help='range of positions for TSS pausing')
    parser.add_argument('--promoter_pause_duration', type=str, default='0,0', help='range for duration of promoter pausing in minutes')
    parser.add_argument("--seed", type=int, help="random seed for reproducibility")

    args = parser.parse_args()

    # Set random seed if specified
    if args.seed is not None:
        random.seed(args.seed)

    # Parse region size or number
    if args.region_size_range and args.region_size_range != 'None':
        min_len, max_len = map(int, args.region_size_range.split(','))
        num_regions = None
    elif args.num_regions:
        min_len, max_len = None, None
        num_regions = args.num_regions
    else:
        raise ValueError("Either --region_size_range or --num_regions must be specified.")

    # Get TSS pausing position range
    promoter_pause_upstream_limit, promoter_pause_downstream_limit = map(int, map(float, args.promoter_pause_position.split(',')))

    # Get TSS pause duration
    promoter_pause_duration_min, promoter_pause_duration_max = map(float, args.promoter_pause_duration.split(','))
    TSS_pause_duration = random.uniform(promoter_pause_duration_min, promoter_pause_duration_max)

    # Get other parameters
    elong_rate_range_min, elong_rate_range_max = map(float, args.elong_rate_range.split(','))
    pause_length_min, pause_length_max = map(float, args.pause_time.split(','))
    pause_duration = random.uniform(pause_length_min, pause_length_max)
    pause_elong_rate = 1 / pause_duration if pause_duration > 0 else float('inf')

    # Read and process input file
    df = pd.read_csv(args.tsv, sep='\t', comment='#')
    df.columns = ['chr', 'start', 'end', 'gene_id', 'feature', 'position', 'strand', 'sequence']

    basename = os.path.basename(args.tsv)
    base_filename = basename.split('.')[0]
    rate_per_gene_dir = os.path.join(args.o, "rate_per_gene")
    os.makedirs(rate_per_gene_dir, exist_ok=True)

    output_filename_regions = os.path.join(rate_per_gene_dir, f"{base_filename}_VariableElongationRateRegions.gtf")

    dna_sequence = getDNAseq(df.sequence)

    df_regions = defineElongRateRegions(base_filename, df, dna_sequence,
                                        min_len, max_len,
                                        elong_rate_range_min, elong_rate_range_max,
                                        pause_elong_rate,
                                        num_regions=num_regions,
                                        num_pauses=args.num_pauses,
                                        promoter_pause_upstream_limit=promoter_pause_upstream_limit,
                                        promoter_pause_downstream_limit=promoter_pause_downstream_limit,
                                        TSS_pause_duration=TSS_pause_duration)

    saveGTF(df_regions, output_filename_regions)

    # Calculate Total Pause Time. Pauses are defined as 1-nucleotide regions.
    pause_regions_df = df_regions[df_regions['region_start_coord'] == df_regions['region_end_coord']]
    if not pause_regions_df.empty:
        total_pause_time_minutes = pause_regions_df['time_to_traverse'].sum()
    else:
        total_pause_time_minutes = 0.0

    output_filename_pausing = os.path.join(rate_per_gene_dir, f"{base_filename}_total_pausing_time.tsv")
    pause_time_df = pd.DataFrame({'total_pause_time_minutes': [total_pause_time_minutes]})
    pause_time_df.to_csv(output_filename_pausing, sep='\t', index=False)
    
    output_filename_nttraversaltime = os.path.join(rate_per_gene_dir, f"{base_filename}_RatesandTraversalTimes.gtf")
    df_nt = ntTraversalTime(df_regions)
    saveGTF(df_nt, output_filename_nttraversaltime)

    print(f"Processing complete. Output files are in: {rate_per_gene_dir}")


