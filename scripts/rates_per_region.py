## Get the region specific elongation rates and the nt level traversal rates
###

import pandas as pd
import random
import argparse
import math
import os

# first set of functions - to be called independently
# define functions:
# 0: function to save gtf files in functions 1, 2, and 5
def saveGTF(input_df, output_filename):
    input_df.to_csv(output_filename, sep = '\t', index = False)

# 1: Create the DNA sequence from the dataframe, concatenating sequence strings for each feature into one string with the entire gene's sequence
# Called like so for each .tsv file: dna_sequence = getDNAseq(df.sequence)
def getDNAseq(dna_sequence):
    return ''.join(str(i) for i in dna_sequence)

# 2: define region-specific elongation rates for each gene, with the goal of making a new .gtf for the gene:
# input is concatenated dna_sequence from gtf file
# called like so for each .gtf file: regions = defineElongRateRegions(input_df_basename, df, dna_sequence)  
# columns of output .gtf will be: chromosome, source, region_start_coord, region_end_coord, rate_initial, rate_final, delta_rate, time_to_traverse, sequence
def defineElongRateRegions(input_df_basename, df, dna_sequence, min_length, max_length, min_rate, max_rate, pause_rate, pause_chance):
    region_data = [] # List to store the start and end coordinates for each region
    gene_length = len(dna_sequence)
    start_coord = 1

    # Make sure region size range is in integer form
    min_length = int(min_length)
    max_length = int(max_length)

    # Make sure elongation rate range is in integer form
    min_rate = int(min_rate)
    max_rate = int(max_rate)
    
    # grab info from the input gtf for the output version
    chromosome = df.iloc[0,0]
    strand = df.iloc[0,6]
    
    abs_start = df.iloc[0, 1]
    abs_end = df.iloc[0, 2]
    
    # define initial rates for the first region, between 0.5 kb/min and 5 kb/min (now in 500 - 5000 nt/min)
    rate_initial = round(random.uniform(min_rate, max_rate), 4)
    
    # parse booleans from input:
    # if true, gene_level will specify that the entire gene will be the "region" (n = 1 region)
    gene_level_rates = args.gene_level
    # if true, flat_rate will make rates invariant over each region
    flat_rates = args.flat_rates

    while start_coord < gene_length: # start going through gene to define regions
        # randomly decide whether to create a pause or a regular region
        if random.random() < pause_chance and start_coord > 1 :
            # create a pol II pause on 1nt
            region_length = 1
            end_coord = start_coord
            sequence = dna_sequence[start_coord - 1 : end_coord]
            
            # make the rate_final very small
            pause_rate_initial = pause_rate
            pause_rate_final = pause_rate
            rate_change_per_nt = (pause_rate_final - pause_rate_initial) / region_length
            
            # add a time_to_traverse for sanity checks
            time_to_traverse = (end_coord - start_coord + 1) / ((rate_final + rate_initial) / 2) 
            
            # add it to the region_data list
            region_data.append((chromosome, abs_start,abs_end, start_coord, end_coord, pause_rate_initial, pause_rate_final, rate_change_per_nt, time_to_traverse, sequence))
            

            # update starting coordinates
            start_coord += region_length
            
        else:                        
            if gene_level_rates == True:
                        region_length = gene_length
            else:
                        # generate a random length for the region
                        region_length = random.randint(min_length, max_length)
            
            # determine end coordinate, making sure it doesn't exceed the total gene sequence length
            end_coord = min((start_coord - 1) + region_length, gene_length)
            
            # extract the sequence for these coordinates
            sequence = dna_sequence[start_coord - 1: end_coord]
            
            if flat_rates == True:
                        rate_final = rate_initial
            else:
                        # determine a random rate_final between 0.5 and 5 kb/min (now in 500 - 5000 nt/min)
                        rate_final = round(random.uniform(min_rate, max_rate), 4)
            
            # determne the change in rate per nt, since it is a linear change (slope) - use the coordinates since the last region was getting messed up by region length
            rate_change_per_nt = (rate_final - rate_initial) / (end_coord - start_coord + 1)
            
            # adjust the rate_final for indexing reasons
            # if you don't do this, the last nt in the region will never reach the rate_final when you increment by rate_change_per_nt
            rate_final = round(rate_final - rate_change_per_nt, 4)
            
            # calculate a time_to_traverse for sanity checks
            time_to_traverse = (end_coord - start_coord + 1) / ((rate_final + rate_initial) / 2) 
            
            # add the start and end coordinates to the region_coordinates list
            #region_data.append((chromosome, source, start_coord, end_coord, rate_initial, rate_final, rate_change_per_nt, time_to_traverse, sequence))
            region_data.append((chromosome, abs_start,abs_end, start_coord, end_coord, strand,rate_initial, rate_final, rate_change_per_nt, time_to_traverse, sequence))
            # update the starting coordiante for the next region
            start_coord = end_coord + 1
            
            if flat_rates == True:
                        rate_initial = round(random.uniform(min_rate, max_rate), 4)
            else:
                        # update the rate_initial for the subseqent region
                        #rate_initial = round(rate_final * random.uniform(0.9,1.1), 3)
                        rate_initial = rate_final
            
    
    # construct the dataframe for the gtf ouptut
    columns = ["chromosome", 'absolute_start','absolute_end', "region_start_coord", "region_end_coord", 'strand',"rate_initial", "rate_final", "rate_change_per_nt", "time_to_traverse", "sequence"]
    regions_df = pd.DataFrame(region_data, columns = columns)
    
    abs_pos = df.iloc[0, 1] if strand == '+' else df.iloc[0, 2]
    
    if strand == '+':
        regions_df["absolute_start"] = regions_df["absolute_start"] - 1 + regions_df["region_start_coord"]
        regions_df["absolute_end"] = regions_df["absolute_start"] - 1 + regions_df["region_end_coord"]
    else:
        regions_df["absolute_start"] = regions_df["absolute_end"] - regions_df["region_end_coord"] + 1
        regions_df["absolute_end"] = regions_df["absolute_end"] - regions_df["region_start_coord"] + 1

    return regions_df

#3: obtain information on the rate and time to traverse for each nt, used for finding stop sites
# called like so for each region_data df: ntTraversalTime(region_data_df)
def ntTraversalTime(region_specific_eRates_df): # base_filename: need to add the file basename here as input for saving file
    
    # grab info from the input gtf for the output version
    chromosome = region_specific_eRates_df.iloc[0,0]
    strand = region_specific_eRates_df.iloc[0,5]
    abs_pos = region_specific_eRates_df.iloc[0, 1] if strand == '+' else region_specific_eRates_df.iloc[0, 2]

    # initialize the array to store the information on how long it takes each nt to be traversed
    traversal_time_per_nt = []
    num_regions = len(region_specific_eRates_df)
    
    # initialize a cumulative nt position for the output df
    cumulative_nt_position = 1
    
    for region in range(num_regions):
        
        # grab what we need from the region specific elongation rate gtf file
        region_length = region_specific_eRates_df['region_end_coord'].iloc[region] - region_specific_eRates_df['region_start_coord'].iloc[region] + 1
        rate_initial = float(region_specific_eRates_df['rate_initial'].iloc[region])
        rate_final = float(region_specific_eRates_df['rate_final'].iloc[region]) 
        rate_change_per_nt = float(region_specific_eRates_df['rate_change_per_nt'].iloc[region])

        
        for nt_position in range(region_length):
            
            rate_for_this_nt = rate_initial + (rate_change_per_nt * nt_position)
            
            time_for_this_nt = 1 / (rate_for_this_nt) # rate is in nt/min
            
            traversal_time_per_nt.append((chromosome, abs_pos, strand,region+1, cumulative_nt_position, time_for_this_nt, rate_for_this_nt, rate_change_per_nt))
            
            # increment the cumulative nt position
            cumulative_nt_position += 1
            
    # construct the dataframe for the gtf ouptut
    columns_nt_table = ["chromosome", "absolute_position", 'strand',"region_number", "nucleotide_coord", "time_for_this_nt", "rate_for_this_nt", "rate_change_per_nt"]
    nt_traversal_time_df = pd.DataFrame(traversal_time_per_nt, columns = columns_nt_table)
    
    if strand == '+':
        nt_traversal_time_df["absolute_position"] = nt_traversal_time_df["absolute_position"] - 1 + nt_traversal_time_df["nucleotide_coord"]
    else:
        nt_traversal_time_df["absolute_position"] = nt_traversal_time_df["absolute_position"] + 1 - nt_traversal_time_df["nucleotide_coord"]
    
    return nt_traversal_time_df

# now we call the main script
if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument('--tsv', default='tmp.tsv', help='full path to input tsv file')
    parser.add_argument('--region_size_range', type=str, default='100,5000', help='range of region sizes (in nt) to be randomly generated', required=False)
    parser.add_argument('--elong_rate_range', type=str, default='500,5000', help='range of elongation rates (in nt/min) to be applied to random regions', required=False)
    parser.add_argument('--pause_occur_probability', type=float, default=0, help='probability of having a pausing event across the isoform', required=False)
    parser.add_argument('--pause_time', default="0.1,0.5", type=str, help="length of the pausing event in minutes", required=False)
    parser.add_argument('--o', type=str, default='./', help='output path')
    parser.add_argument('--flat_rates', action='store_true', help='all nucleotides within the region have the same elongation rate. If used with gene_level all nucleotides will have the same elongation rates.')
    parser.add_argument('--gene_level', action='store_true', help='get only one region across the gene. If used with flat_rates all nucleotides will have the same elongation rates.')
    args = parser.parse_args()

    # set variables for arguments needed for functions 1, 2, 3
    #region_size_range = [int(x) for x in args.region_size_range.split(',')]
    region_size_range_min, region_size_range_max = map(int, args.region_size_range.split(','))

    #elong_rate_range = [int(x) for x in args.elong_rate_range.split(',')]
    elong_rate_range_min, elong_rate_range_max = map(float, args.elong_rate_range.split(','))

    #Get pause duration and convert it to elongation rate for calculation
    pause_chance = args.pause_occur_probability
    pause_length_min, pause_length_max = map(float, args.pause_time.split(','))
    pause_duration = random.uniform(pause_length_min, pause_length_max)
    pause_elong_rate=1/pause_duration

    
    # create strings for arguments needed in functions 1, 2, and 3 for output file naming
    region_size_name = str(region_size_range_min) + '-' + str(region_size_range_max)
    elong_rate_name = str(elong_rate_range_min) + '-' + str(elong_rate_range_max)
    pause_rate_name = str(pause_elong_rate)
    pause_chance_name = str(pause_chance)
    
    # read in the input .tsv file
    df = pd.read_csv(args.tsv, sep='\t', comment='#')
    df.columns = ['chr','start','end','gene_id','feature','position','strand','sequence']

    # Generate output name and file for output file with read and substitution info
    base_filename = os.path.splitext(os.path.splitext(os.path.basename(args.tsv))[0])[0]

    #Generate dir to store output
    rate_per_gene_dir = os.path.join(args.o, "rate_per_gene")
    os.makedirs(rate_per_gene_dir, exist_ok=True)


    # name the region specific elongation rate .gtf: concatenate the base filename, the region size range, elongation rate range (nt/min), pause rate (nt/min), and pause chance per random region generated
    output_filename_regions = str(rate_per_gene_dir) + '/'+base_filename + '_VariableElongationRateRegions.gtf'

    # use fxn #1 to get the string of the gene DNA sequence
    dna_sequence = getDNAseq(df.sequence)

    # use fxn #2 to make the region-specific elongation rate file
    df_regions = defineElongRateRegions(base_filename, df, dna_sequence, region_size_range_min, region_size_range_max, elong_rate_range_min, elong_rate_range_max, pause_elong_rate, pause_chance)
    df.loc['region_number'] = range(1, len(df.columns) + 1)

    # use fxn #0 to save the output gtf for function #2
    saveGTF(df_regions, output_filename_regions)
    
    # name the nt traversal time gtf: use info from output filename regions but change the descriptor
    output_filename_nttraversaltime = str(rate_per_gene_dir) + '/' + base_filename + '_RatesandTraversalTimes.gtf'
    
    # use fxn #3 to make the nt traversal time file and save it
    df_nt = ntTraversalTime(df_regions)
    
    # use fxn #0 to save the output gtf for function #3
    saveGTF(df_nt, output_filename_nttraversaltime)




