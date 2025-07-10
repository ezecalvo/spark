# second set of functions - to be called after the first set

import pandas as pd
import random
import argparse
import math
import os
import string

# 1: Create the DNA sequence from the dataframe, concatenating sequence strings for each feature into one string with the entire gene's sequence
# Called like so for each .gtf file: dna_sequence = getDNAseq(df.sequence)
def getDNAseq(dna_sequence):
    return ''.join(str(i) for i in dna_sequence)

# 4: find nt stop sites for each labeling period
# called like so for each labeling period: findStopnt(df_nt, labeling_duration, startsite), where df_nt is the output of ntTraversalTime

# find the nt position reached up to the labeling duration, given a nt startsite
def findStopnt(nt_traversal_time_df, labeling_duration, startsite):
    
    #initialize cumulative time object
    cumulative_time_elapsed = 0
    
    # initialize variables to store the penultimate loop values for when labeling duration is exceeded
    prev_nt_coord = None
    prev_cumulative_time_elapsed = 0
    
    # grab the nt at which this function will begin
    startsite_nt_coord = nt_traversal_time_df['nucleotide_coord'].iloc[startsite - 1] #startsite was in 1 index, shift back to 0 based for this
    
    # find the point at which 5 minutes will have elapsed (but not exceed labeling duration)
    for nt in range(startsite_nt_coord - 1, len(nt_traversal_time_df)):
        
        # start by storing the current nucleotide coordinate and cumulative time before adding the new one
        prev_nt_coord = nt + 1 # for 1-based indexing
        prev_cumulative_time_elapsed = cumulative_time_elapsed 

        # add the time per nt value to cumulative time elapsed for this nucleotide
        cumulative_time_elapsed += nt_traversal_time_df.loc[nt, 'time_for_this_nt']
        
        # if cumulative time elapsed exceeds the labeling period, stop the loop and backtrack one nt 
        if cumulative_time_elapsed > labeling_duration:
            endsite_nt_coord = prev_nt_coord
            cumulative_time_elapsed = prev_cumulative_time_elapsed
            break
        
        # save the coordinate of the nt at which the cumulative time elapsed exceeded the labeling duration
        endsite_nt_coord = nt + 1
        
        # save the rate at that nucleotide for setting up the next start site
        #endsite_nt_rate = nt_traversal_time_df['rate_for_this_nt'].iloc[nt + 1] # for 1-based indexing
    
    #return(startsite_nt_coord, endsite_nt_coord, endsite_nt_rate)
    return(endsite_nt_coord)


# 5: Perform nucleotide conversions within the labeled region
def convert_base(read, nt_incorporation_rate, nt_conversion_rate, start_site, stop_site, from_base, to_base):
    
    # this will read the sequence of the region created during the labeling period in question (5eU, 6sG, or 4sU) and put it into a substring
    labeled_region = read[start_site:stop_site]
    
    # get all of the positions in the substring that correspond to the "from_base" character
    pos = [pos for pos, char in enumerate(labeled_region) if char == from_base]
    
    # determine how many of the from_base bases will have modified nts incorporated 
    num_bases_incorporated = int(math.ceil(len(pos) * nt_incorporation_rate))
    
    # collect all the positions that will have modified nucleotides incorporated
    pos_conv = random.sample(pos, num_bases_incorporated)
    
    # convert string of the nt read to a list so nt characters can be changed
    conv_read = list(read)
    
    # actually substitute the nts in the read using this loop
    for i in pos_conv:
        
        # this needs to take into account the actual conversion rate for the given nucleotide (N/A for 5eU, 0.85-0.90 for 4sU, 0.35-0.55 for 6sG)
        # it's definitely going to slow down everything though, so maybe it's not the right way to implement it...
        if random.random() < random.uniform(nt_conversion_rate[0], nt_conversion_rate[1]):
            
            conv_read[start_site + i] = to_base
            
    # as curretnly written, this will output the positions where the nt was incorporated, but not necessarily converted. Do I need to update this to also output the positions with real conversions?
    return ''.join(conv_read), nt_incorporation_rate, [start_site + i for i in pos_conv]


# 6: Mutate random positions in the mRNA to simulate sequencing errors
def add_seq_errs(read, seq_error_rate_range):
    
    # get the length of the total transcript made across all 3 labeling periods
    read_length = len(read)
    
    # determine the sequencing error rate for this transcript
    seq_err_rate = round(random.uniform(seq_error_rate_range[0], seq_error_rate_range[1]), 5)
    
    # figure out how many nts should be sequenced in error
    num_seq_err_pos = math.ceil(read_length * seq_err_rate)
    
    # index the read
    pos = [pos for pos, char in enumerate(read)]
    
    # randomly select positions to sequence in error
    pos_seq_err = random.sample(pos, int(num_seq_err_pos))
    
    # sort the positions in ascending order to make sure the errors are applied left-to-right in the sequence
    pos_seq_err.sort()
    
    # covert the read to a list that can have nts changed
    seq_err_read = list(read)
    
    #iterate through the positions to be sequenced in error
    for i in pos_seq_err:
        
        # grab a random nt
        rand_nt = random.choice(["A", "C", "G", "T"])
        
        while (rand_nt == seq_err_read[i]):
            
            # if the random nt is the same as the one at that position, grab a new one. rinse and repeat until they are different
            rand_nt = random.choice(["A", "C", "G", "T"]) # weight probability of changing T to C, etc, but keep selection of nt to be changed as random
        
        # I'm not sure what this is for...
        ref_nt = seq_err_read[i]
        
        # again, not sure why this is necessary
        seq_err_nt = rand_nt
        
        # swap the nt at this position with the incorrectly sequenced nt
        seq_err_read[i] = seq_err_nt
        
    # Convert the list back into a string
    seq_err_read = "".join(seq_err_read)

    # return the read with errors, the rate of errors, and the positions where sequencing errors occurred
    return seq_err_read, seq_err_rate, pos_seq_err


# 7: lastly, create the output name for the gtf to feed into the chopper script
def createOutputName(gtf_file, label_time, bkgd_perc_range, u2c_perc_range, num_reads):
    
    gtf_name = gtf_file.strip(".gtf")
    
    output_file_name = str(gtf_name) + "_l" + str(label_time) + "_b" + str(bkgd_perc_range) + "_p" + str(u2c_perc_range) + "_reads" + str(num_reads) + ".csv.gz"

#8: add molecule id (this will also be the read name for long reads)
def id_generator(size=12, chars=string.ascii_uppercase + string.digits):
    return "@" + ''.join(random.choice(chars) for _ in range(size - 1))   

# now we call the main script
if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument('--region_file', default='tmp.gtf', help='Input region specific elongation rate GTF file name')
    parser.add_argument('--nt_file', default='tmp.gtf', help='Input ntTraversalTime GTF file name')
    parser.add_argument('--l', type=int, default=5, help='Labeling time in minutes')
    parser.add_argument('--seq_err', default='0.0001,0.0002', help='Illumina sequencing error rate range')
    parser.add_argument('--nt_inc_rate', type=float, default=0.1, help = 'modified nucleotide INCORPORATION rate for 5eU, 4sU, and 6sG', required=False)
    parser.add_argument('--subs_5eU', type=str, default='1,1', help='range of percent CONVERSION for 5eU', required=False)
    parser.add_argument('--subs_6sG', type=str, default='0.35,0.55', help='range of percent CONVERSION for 6sG', required=False)
    parser.add_argument('--subs_4sU', type=str, default='0.85,0.90', help='range of percent CONVERSION for 4sU', required=False)
    parser.add_argument('--reads', type=int, default=200, help='Number of reads to generate')
    parser.add_argument('--o', type=str, default='./', help='output path')
    parser.add_argument('--treatment', type=str, default='no DRB', help='usage of DRB or not')
    args = parser.parse_args()

    # set variables for arguments needed for functions 4-7
    nt_incorporation_rate = args.nt_inc_rate
    u2u_percent_range = list(map(lambda x: float(x), args.subs_5eU.split(',')))
    g2a_percent_range = list(map(lambda x: float(x), args.subs_6sG.split(',')))  # *** ORIGINAL ***
    u2c_percent_range = list(map(lambda x: float(x), args.subs_4sU.split(',')))    
    seqerr_range = list(map(lambda x: float(x), args.seq_err.split(',')))
    
    # create strings for arguments needed in functions 4-7 for output file naming
    nt_inc_rate_name = str(int(nt_incorporation_rate*100)) + '-' + str(int(nt_incorporation_rate*100)) # convert to percent for compatibility with later scripts
    u2u_subrate_name = str(u2u_percent_range[0]) + '-' + str(u2u_percent_range[1])
    g2a_subrate_name = str(g2a_percent_range[0]) + '-' + str(g2a_percent_range[1])
    u2c_subrate_name = str(u2c_percent_range[0]) + '-' + str(u2c_percent_range[1])
    seqerr_name = str((seqerr_range[0]*100)) + '-' + str((seqerr_range[1]*100)) # convert to percent without decimals

    # read in the input region specific elongation rate .gtf file
    df = pd.read_csv(args.region_file, sep='\t', comment='#')
    df.columns = ["chromosome", 'absolute_start','absolute_end', "region_start_coord", "region_end_coord", 'strand',"rate_initial", "rate_final", "rate_change_per_nt", "time_to_traverse", "sequence"]
    
    # read in the input nt traversal time .gtf file
    df_nt = pd.read_csv(args.nt_file, sep='\t', comment='#')
    df_nt.columns = ["chromosome", "absolute_position", 'strand',"region_number", "nucleotide_coord", "time_for_this_nt", "rate_for_this_nt", "rate_change_per_nt"]
    
    # use fxn #1 to get the string of the gene DNA sequence
    dna_sequence = getDNAseq(df.sequence)

    # Generate output name and file for output file with read and substitution info
    # Generate output name and file for output file with read and substitution info
    filename = os.path.splitext(os.path.basename(args.region_file))[0]
    # Extract only the gene ID (first part before the first underscore)
    base_filename = filename.split("_")[0]
    # Construct the output filename
    output_filename = str(args.o) + base_filename + "_KB.tsv"

    # open the output file for editing
    output_file = open(output_filename, 'w')
    
    #Generate initiation times for generating reads, one for each read spanning the entire labeling experiment (labeling time x 3)
    # initiation time is treated as the time remaining in the experiment, for example if initiation = 2, then there is 2 minutes remaining in the experiment
    initiation_times=[]
    for i in range(0,args.reads):
        initiation_temp=random.uniform(0,args.l*3) 
        initiation_times.append(initiation_temp)
    
    # initialize the read number counter
    n=0

    # initialize the master data list that will be output at the end
    data_list = []
    max_transcript_length = df["region_end_coord"].max()
    strand=df_nt.loc[0, 'strand']

    # the underscore here means the value has no import on the function (here it is the read number)
    for _ in range(args.reads):
        # grab the nth initiation time
        test = initiation_times[n]
        initiation = test
        
        # set the startsite for 5eU labeling (label #1)
        starts_5eU_pos = 1 if args.treatment == "DRB" else random.randint(1, max_transcript_length)       
        
        # findStopnt inputs: findStopnt(nt_traversal_time_df, labeling_duration, startsite):
        # findStopnt outputs: return(endsite_nt_coord)
        
        # use function #4 to find the stop site nt coordinate for 5eU labeling (the first 5 minutes)
        # use initiation as the labeling duration, as 5eU will be available for the entirety of the experiment (if it doesn't reach the end of the gene, then that transcript wouldn't be fully transcribed)
        stop_site_5eU = findStopnt(df_nt, initiation, starts_5eU_pos) 
        
        if initiation > args.l*2: # if there are more than 2 timepoints worth of minutes remaining
            # with more than 2 timepoints worth of minutes remaining, find the position where 6sG will begin labeling
            starts_6sG_pos = findStopnt(df_nt, initiation-(args.l*2), 1)
            # using the 6sG start site, find the stop site for 6sG, given that it will be available for the final 2 labeling timepoints
            stop_site_6sG = findStopnt(df_nt, args.l*2, starts_6sG_pos)
            
        else: # meaning if initiation is not greater than 10, 6sG will be available for the duration of that mRNA's transcription
            # transcription is beginning for this mRNA while both 5eU and 6sG are available
            starts_6sG_pos = 1
            # find the stop site for 6sG from the beginning of the gene, given that initiation is all the remaining time in the experiment (if it doesn't reach the end of the gene, then that transcript wouldn't be fully transcribed)
            stop_site_6sG = findStopnt(df_nt, initiation, starts_6sG_pos)

        if initiation>args.l: # if there is more than 1 timepoints worth of minutes remaining
            # with more than 1 timepoints worth of minutes remaining, find the position where 4sU will begin labeling
            starts_4sU_pos = findStopnt(df_nt, initiation-args.l, 1)
            # using the 4sU start site, find the stop site for 4sU, given that it will be available for the final labeling timepoint
            stop_site_4sU = findStopnt(df_nt, args.l, starts_4sU_pos)#4sU will be available for the entire 3 timepoint
            
        else: # meaning if initiation is not greater than 5, 4sU will be available for the duration of that mRNA's transcription
            # transcription is beginning for this mRNA while 5eU, 6sG, and 4sU are all available
            starts_4sU_pos = 1
            # find the stop site for 4sU from the beginning of the gene, given that initiation is all the remaining time in the experiment (if it doesn't reach the end of the gene, then that transcript wouldn't be fully transcribed)
            stop_site_4sU = findStopnt(df_nt, initiation, starts_4sU_pos)

        # construct the read sequence, using the coordinates identified from 5eU labeling, which was available for the entirety of transcription
        read = dna_sequence[starts_5eU_pos:int(stop_site_5eU)]
        
        # convert 5eU positions
        converted_read, percent_u2u, converted_positions_5eU = convert_base(read, nt_incorporation_rate, u2u_percent_range, int(starts_5eU_pos), int(stop_site_5eU),'T','T')
        
        # convert 6sG positions
        converted_read, percent_g2a, converted_positions_6sG = convert_base(converted_read, nt_incorporation_rate, g2a_percent_range, int(starts_6sG_pos), int(stop_site_6sG),'G','A')
        
        # convert 4sU positions
        converted_read, percent_t2c, converted_positions_4sU = convert_base(converted_read, nt_incorporation_rate, u2c_percent_range , int(starts_4sU_pos), int(stop_site_4sU),'T','C')

        #bg_mutated_read, percentage_bg_mutations, mutated_bg_positions = mutate_bg(converted_read, list(map(float, args.b.split(','))), int(starts_5eU_pos), int(stop_site_5eU))
        seq_err_read, percentage_seq_err, seq_err_positions = add_seq_errs(converted_read, list(map(float, args.seq_err.split(','))))

        molecule_id=id_generator()

        n = n + 1
        row_data = [initiation, molecule_id,strand,percent_u2u, percent_g2a, percent_t2c, starts_5eU_pos, stop_site_5eU, starts_6sG_pos, stop_site_6sG, starts_4sU_pos, stop_site_4sU, percentage_seq_err, converted_positions_5eU, converted_positions_6sG, converted_positions_4sU, seq_err_positions, seq_err_read]
      #  row_data.columns = ["initiation_time", 'percent_u2u', 'percent_g2a', 'percent_t2c', 'starts_5eU_pos', 'stop_site_5eU', 'starts_6sG_pos', 'stop_site_6sG', 'starts_4sU_pos', 'stop_site_4sU', 'percentage_seq_err', 'converted_positions_5eU', 'converted_positions_6sG', 'converted_positions_4sU', 'seq_err_positions','full_molecule_sequence']

        # Add the row_data to the data_list
        data_list.append(row_data)

        # Convert to DataFrame with column names

    df_for_export = pd.DataFrame(data_list, columns=["initiation_time", 'molecule_id','strand','percent_u2u', 'percent_g2a', 'percent_t2c', 'starts_5eU_pos',
        'stop_site_5eU', 'starts_6sG_pos', 'stop_site_6sG', 'starts_4sU_pos', 'stop_site_4sU', 'percentage_seq_err',
        'converted_positions_5eU', 'converted_positions_6sG', 'converted_positions_4sU', 'seq_err_positions','full_molecule_sequence'])

    # Save to a tab-separated file
    df_for_export.to_csv(output_filename, sep='\t', index=False)


