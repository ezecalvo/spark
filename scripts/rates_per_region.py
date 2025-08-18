import pandas as pd
import random
import argparse
import math
import os

def saveGTF(input_df, output_filename):
    input_df.to_csv(output_filename, sep='\t', index=False)

def getDNAseq(dna_sequence):
    return ''.join(str(i) for i in dna_sequence)

def defineElongRateRegions(input_df_basename, df, dna_sequence,
                           min_length, max_length,
                           min_rate, max_rate,
                           pause_rate, pause_chance,
                           num_regions=None):

    region_data = []
    gene_length = len(dna_sequence)
    start_coord = 1

    min_rate = int(min_rate)
    max_rate = int(max_rate)

    chromosome = df.iloc[0, 0]
    strand = df.iloc[0, 6]
    abs_start = df.iloc[0, 1]
    abs_end = df.iloc[0, 2]

    rate_initial = round(random.uniform(min_rate, max_rate), 4)

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

    for region_length in lengths:
        end_coord = start_coord + region_length - 1
        sequence = dna_sequence[start_coord - 1:end_coord]

        if random.random() < pause_chance and start_coord > 1 and region_length > 1:
            pause_rate_initial = pause_rate
            pause_rate_final = pause_rate
            rate_change_per_nt = 0
            time_to_traverse = region_length / pause_rate

            region_data.append((chromosome, abs_start, abs_end, start_coord, end_coord, strand,
                                pause_rate_initial, pause_rate_final, rate_change_per_nt, time_to_traverse, sequence))
        else:
            flat_rates = args.flat_rates
            if flat_rates:
                rate_final = rate_initial
            else:
                rate_final = round(random.uniform(min_rate, max_rate), 4)

            rate_change_per_nt = (rate_final - rate_initial) / region_length
            rate_final = round(rate_final - rate_change_per_nt, 4)
            time_to_traverse = region_length / ((rate_final + rate_initial) / 2)

            region_data.append((chromosome, abs_start, abs_end, start_coord, end_coord, strand,
                                rate_initial, rate_final, rate_change_per_nt, time_to_traverse, sequence))

            rate_initial = rate_final if not flat_rates else round(random.uniform(min_rate, max_rate), 4)

        start_coord = end_coord + 1

    columns = ["chromosome", 'absolute_start', 'absolute_end', "region_start_coord", "region_end_coord",
               'strand', "rate_initial", "rate_final", "rate_change_per_nt", "time_to_traverse", "sequence"]
    regions_df = pd.DataFrame(region_data, columns=columns)

    if strand == '+':
        regions_df["absolute_start"] = regions_df["absolute_start"] - 1 + regions_df["region_start_coord"]
        regions_df["absolute_end"] = regions_df["absolute_start"] - 1 + regions_df["region_end_coord"]
    else:
        regions_df["absolute_start"] = regions_df["absolute_end"] - regions_df["region_end_coord"] + 1
        regions_df["absolute_end"] = regions_df["absolute_end"] - regions_df["region_start_coord"] + 1

    return regions_df

def ntTraversalTime(region_specific_eRates_df):
    chromosome = region_specific_eRates_df.iloc[0, 0]
    strand = region_specific_eRates_df.iloc[0, 5]
    abs_pos = region_specific_eRates_df.iloc[0, 1] if strand == '+' else region_specific_eRates_df.iloc[0, 2]

    traversal_time_per_nt = []
    cumulative_nt_position = 1

    for region in range(len(region_specific_eRates_df)):
        region_length = (region_specific_eRates_df['region_end_coord'].iloc[region] -
                         region_specific_eRates_df['region_start_coord'].iloc[region] + 1)
        rate_initial = float(region_specific_eRates_df['rate_initial'].iloc[region])
        rate_change_per_nt = float(region_specific_eRates_df['rate_change_per_nt'].iloc[region])

        for nt_position in range(region_length):
            rate_for_this_nt = rate_initial + (rate_change_per_nt * nt_position)
            time_for_this_nt = 1 / rate_for_this_nt

            traversal_time_per_nt.append((chromosome, abs_pos, strand, region + 1,
                                           cumulative_nt_position, time_for_this_nt,
                                           rate_for_this_nt, rate_change_per_nt))
            cumulative_nt_position += 1

    columns_nt_table = ["chromosome", "absolute_position", 'strand', "region_number", "nucleotide_coord",
                        "time_for_this_nt", "rate_for_this_nt", "rate_change_per_nt"]
    nt_traversal_time_df = pd.DataFrame(traversal_time_per_nt, columns=columns_nt_table)

    if strand == '+':
        nt_traversal_time_df["absolute_position"] = (nt_traversal_time_df["absolute_position"] - 1 +
                                                       nt_traversal_time_df["nucleotide_coord"])
    else:
        nt_traversal_time_df["absolute_position"] = (nt_traversal_time_df["absolute_position"] + 1 -
                                                       nt_traversal_time_df["nucleotide_coord"])

    return nt_traversal_time_df

if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument('--region_size_range', type=str, help='min,max range of region sizes (in nt)')
    parser.add_argument('--num_regions', type=int, help='desired number of regions to split the gene into')

    parser.add_argument('--elong_rate_range', type=str, default='500,5000', help='range of elongation rates (nt/min)')
    parser.add_argument('--pause_occur_probability', type=float, default=0, help='probability of pausing event')
    parser.add_argument('--pause_time', type=str, default="0.1,0.5", help='length of pausing in minutes')
    parser.add_argument('--tsv', type=str, default='tmp.tsv', help='input TSV file')
    parser.add_argument('--o', type=str, default='./', help='output path')
    parser.add_argument('--flat_rates', action='store_true', help='make rates constant within regions')
    parser.add_argument('--gene_level', action='store_true', help='use one region per gene')

    args = parser.parse_args()

    #Parse region size or number
    # Convert 'None' string to real None
    if args.region_size_range == 'None':
        region_size_range = None
    if args.num_regions == 'None':
        num_regions = None

    if args.region_size_range:
        min_len, max_len = map(int, args.region_size_range.split(','))
        num_regions = None
    else:
        min_len, max_len = None, None
        num_regions = args.num_regions


    elong_rate_range_min, elong_rate_range_max = map(float, args.elong_rate_range.split(','))
    pause_length_min, pause_length_max = map(float, args.pause_time.split(','))
    pause_duration = random.uniform(pause_length_min, pause_length_max)
    pause_elong_rate = 1 / pause_duration

    df = pd.read_csv(args.tsv, sep='\t', comment='#')
    df.columns = ['chr', 'start', 'end', 'gene_id', 'feature', 'position', 'strand', 'sequence']

    base_filename = os.path.splitext(os.path.splitext(os.path.basename(args.tsv))[0])[0]
    rate_per_gene_dir = os.path.join(args.o, "rate_per_gene")
    os.makedirs(rate_per_gene_dir, exist_ok=True)

    output_filename_regions = os.path.join(rate_per_gene_dir, f"{base_filename}_VariableElongationRateRegions.gtf")

    dna_sequence = getDNAseq(df.sequence)

    df_regions = defineElongRateRegions(base_filename, df, dna_sequence,
                                        min_len, max_len,
                                        elong_rate_range_min, elong_rate_range_max,
                                        pause_elong_rate, args.pause_occur_probability,
                                        num_regions=num_regions)

    saveGTF(df_regions, output_filename_regions)

    output_filename_nttraversaltime = os.path.join(rate_per_gene_dir, f"{base_filename}_RatesandTraversalTimes.gtf")
    df_nt = ntTraversalTime(df_regions)
    saveGTF(df_nt, output_filename_nttraversaltime)

