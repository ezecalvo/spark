import pandas as pd
import random
import argparse
import math
import os
import string
import re
import numpy as np

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
        original_time = nt_traversal_time_df.loc[nt, 'time_for_this_nt']
        jittered_time = original_time * random.uniform(0.9, 1.1)
        cumulative_time_elapsed += jittered_time
        
        # if cumulative time elapsed exceeds the labeling period, stop the loop and backtrack one nt 
        if cumulative_time_elapsed > labeling_duration:
            endsite_nt_coord = prev_nt_coord
            cumulative_time_elapsed = prev_cumulative_time_elapsed
            break
        
        # save the coordinate of the nt at which the cumulative time elapsed exceeded the labeling duration
        endsite_nt_coord = nt + 1
        
    return(endsite_nt_coord)


# 5: Perform nucleotide conversions within the labeled region
def convert_base(read, nt_incorporation_rate, nt_conversion_rate, start_site, stop_site, from_base, to_base=None):
    # Extract the labeled region of the read
    labeled_region = read[start_site:stop_site]
    
    # Find positions of all 'from_base' in the labeled region
    pos = [pos for pos, char in enumerate(labeled_region) if char == from_base]
    
    nt_incorporation_rate_to_use = random.uniform(nt_incorporation_rate[0], nt_incorporation_rate[1])
    # Determine how many will be incorporated
    num_bases_incorporated = int(math.ceil(len(pos) * nt_incorporation_rate_to_use))
    
    # Randomly choose which ones are incorporated
    pos_inc = random.sample(pos, num_bases_incorporated) if num_bases_incorporated > 0 else []
    
    # If to_base is not specified, skip conversion
    if to_base is None or to_base == "":
        num_bases_converted = 0
        pos_to_convert = []
    else:
        # Determine what percent of those will actually be converted
        conversion_fraction = random.uniform(nt_conversion_rate[0], nt_conversion_rate[1])
        num_bases_converted = int(math.floor(len(pos_inc) * conversion_fraction))
        # Choose randomly which of the incorporated bases are converted
        pos_to_convert = random.sample(pos_inc, num_bases_converted) if num_bases_converted > 0 else []
    
    # Convert to list and apply the conversion
    conv_read = list(read)
    for i in pos_to_convert:
        conv_read[start_site + i] = to_base  # start_site is added because `pos` is relative to labeled_region
    
    return ''.join(conv_read), nt_incorporation_rate_to_use, [start_site + i for i in pos_to_convert], [start_site + i for i in pos_inc]




# 6: Mutate random positions in the mRNA to simulate sequencing errors
def add_seq_errs(read, seq_error_rate_range):
    
    # get the length of the total transcript made across all 3 labeling periods
    read_length = len(read)
    
    # determine the sequencing error rate for this transcript
    seq_err_rate = round(random.uniform(seq_error_rate_range[0], seq_error_rate_range[1]), 5)
    
    # figure out how many nts should be sequenced in error
    num_seq_err_pos = np.random.binomial(read_length, seq_err_rate)
    
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


# 7: create the output name for the gtf to feed into the chopper script
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
    parser.add_argument('--experiment_time', type=int, default=15, help='Labeling time in minutes')
    parser.add_argument('--mRNAcount', type=int, default=5000, help='Number of mRNA molecules to simulate')  
    parser.add_argument('--seq_err', default='0.0001,0.0002', help='Illumina sequencing error rate range')
    parser.add_argument('--nt_inc_rate',type=str,default='0.09,0.1',help='Comma-separated range of nucleotide incorporation proportions')
    parser.add_argument('--subs_rate',type=str,default='0.95,1',help='Comma-separated range of nucleotide substitution proportions')
    parser.add_argument('--labeling_base', type=str, default='T', help='Specify the base analog that will get incorporated into the nascentRNA', required=False)
    parser.add_argument('--sub_base', type=str, help='Specify the identity of the base after conversion', required=False)
    parser.add_argument('--initi_rate',type=str,default='10,60',help='Comma-separated range of seconds for initiation rate')
    parser.add_argument('--bkg_molecules', type=float,default=0,help='Proportion of molecules that are background (0 or 1)')
    parser.add_argument('--path_to_tpm', type=str,help='Full path to the TPM per gene file')
    parser.add_argument('--o', type=str, default='./', help='output path')
    parser.add_argument("--drb", action="store_true", help="DRB treatment experiment for transcription synchronization")
    parser.add_argument("--seed", type=int, help="random seed for reproducibility")

    args = parser.parse_args()



    # Set random seed if specified
    if args.seed is not None:
        random.seed(args.seed)

    nt_incorporation_rate = list(map(lambda x: float(x), args.nt_inc_rate.split(',')))
    sub_rate_percent_range = list(map(lambda x: float(x), args.subs_rate.split(',')))
    

    seqerr_range = list(map(lambda x: float(x), args.seq_err.split(',')))

    # read in the input region specific elongation rate .gtf file
    df = pd.read_csv(args.region_file, sep='\t', comment='#')
    df.columns = ["chromosome", 'absolute_start','absolute_end', "region_start_coord", "region_end_coord", 'strand',"rate_initial", "rate_final", "rate_change_per_nt", "time_to_traverse", "sequence"]
    
    # read in the input nt traversal time .gtf file
    df_nt = pd.read_csv(args.nt_file, sep='\t', comment='#')
    df_nt.columns = ["chromosome", "absolute_position", 'strand',"region_number", "nucleotide_coord", "time_for_this_nt", "rate_for_this_nt", "rate_change_per_nt"]
    
    # use fxn #1 to get the string of the gene DNA sequence

    dna_sequence = getDNAseq(df.sequence)

    # Generate output name and file for output file with read and substitution info
    filename = os.path.splitext(os.path.basename(args.region_file))[0]
    # Extract only the gene ID (first part before the first underscore)
    base_filename = filename.split("_")[0]
    # Construct the output filename
    output_filename = str(args.o) + base_filename + ".tsv.gz"

    # open the output file for editing
    output_file = open(output_filename, 'w')

    # initialize the read number counter
    n=0

    # initialize the master data list that will be output at the end
    data_list = []

    #Get expected subs
    from_base=args.labeling_base
    to_base=args.sub_base


    max_transcript_length = df["region_end_coord"].max()

    # Simulate mRNAcount mRNA molecules
    molecule_count = args.mRNAcount

    #For each cell, generate a series of exponentially distributed waiting times (with mean 10 seconds), and accumulate initiation times up to the experiment length
    low, high = map(int, args.initi_rate.split(","))
    initiation_interval_sec = random.randint(low, high)
    initiation_rate_per_sec = 1.0 / initiation_interval_sec
    
    # Calculate how many molecules were already initiated at any given point in time. That is equal to the initiation rate times how long it takes to transverse the gene (length/elong_rate + total paused time)
    #Get total pausing time
    input_path = args.region_file
    pausing_file_path = input_path.replace('_VariableElongationRateRegions.gtf', '_total_pausing_time.tsv')

    total_pausing_times = 0.0 # Default value if file is not found or is empty
    pausing_df = pd.read_csv(pausing_file_path, sep='\t')
    if not pausing_df.empty:
        total_pausing_times = pausing_df['total_pause_time_minutes'].iloc[0]
    #max_transcript_length #This is in nt
    mean_elong_rate=df_nt['rate_for_this_nt'].mean() #This is in nt per min
    n_pol=(initiation_interval_sec/60)*(max_transcript_length/mean_elong_rate+(total_pausing_times))

    #This will make DRB simulate everything
    if args.drb:
        n_pol=0

    #Now I'll scale how many molecules need to initiate
    scaling_factor=molecule_count/(n_pol+(initiation_interval_sec/60)*args.experiment_time)

    already_initiated=n_pol*scaling_factor
    already_initiated=int(already_initiated)
    n_mol_to_initiate=scaling_factor*args.experiment_time*initiation_interval_sec/60

    #This chunk should distribute already_initiated across the gene weighting for the elongation rate of each region
    #For slower regions, we are more likely to find an elongating Pol II
    #inverse_weights = 1 / df_nt['rate_for_this_nt']
    # Use these new inverse weights for sampling
    #chosen_samples = df_nt.sample(n=already_initiated, weights=inverse_weights, replace=True, random_state=42)



    # Identify single-nt regions
    single_nt_mask = df['region_start_coord'] == df['region_end_coord']
    single_nt_df = df[single_nt_mask]

    # Calculate expected count for each based on weights
    inverse_weights = df['time_to_traverse']
    total_inverse_weight = inverse_weights.sum()
    expected_counts = (inverse_weights * already_initiated / total_inverse_weight).round().astype(int)

    # For single_nt regions: use min(scaling_factor, expected_count)
    single_nt_sample_counts = expected_counts[single_nt_mask].clip(upper=int(scaling_factor))

    # Repeat each single-nt region row the appropriate number of times
    single_nt_samples = pd.DataFrame([
        row for idx, count in single_nt_sample_counts.items() for row in [df.loc[idx]] * count
    ])


    # For remaining needed samples, do weighted sampling from the other regions
    n_single_nt = single_nt_sample_counts.sum()
    other_regions = df[~single_nt_mask]
    n_to_sample = already_initiated - n_single_nt
    other_samples = other_regions.sample(n=n_to_sample,weights=inverse_weights[~single_nt_mask],replace=True,random_state=42)
    

    # Combine
    chosen_samples = pd.concat([single_nt_samples, other_samples])

#After weighting, get a random distribution of Pols within each region:
    random_positions = np.random.randint(low=chosen_samples['region_start_coord'], high=chosen_samples['region_end_coord'] + 1)

    start_label_pos_list = random_positions.tolist()    


    #Assign the starting positions
    #start_label_pos_list = []
    #start_label_pos_list=chosen_samples['nucleotide_coord'].tolist()

    initiation_times = []
    experiment_time_sec = args.experiment_time * 60 #Convert experiment time to secs



    #For the molecules at the pause site, give them a bit more time to elongate, this is just to avoid syncronization at the TSS and all of them having to wait the full pausing time
    max_value = int(total_pausing_times * 60)
    for i in range(0, len(single_nt_samples)):
        initiation_temp = experiment_time_sec+random.randint(0, max_value)
        initiation_times.append(initiation_temp)
    

    for i in range(0, len(other_samples)):
        initiation_temp = experiment_time_sec
        initiation_times.append(initiation_temp)



   
    #Now create new initiation events
    # Calculate the average number of initiation events we expect per cell during the experiment
    avg_initiations_per_cell = experiment_time_sec * initiation_rate_per_sec
    #The scaling factor is equal to the number of alleles needed for the desired number of mRNAs
    num_cells_needed=scaling_factor

    # Run the simulation for the calculated number of cells to generate a pool of molecules
    all_new_initiations = []
    if num_cells_needed > 0 and n_mol_to_initiate > 0:
        for cell in range(int(num_cells_needed)):
            initiation_times_cell = []
            current_time = 0.0
            while current_time < experiment_time_sec:
                wait = np.random.exponential(1.0 / initiation_rate_per_sec)
                current_time += wait
                if current_time < experiment_time_sec:
                    initiation_times_cell.append(current_time)
            all_new_initiations.extend(initiation_times_cell)

    final_new_initiations = all_new_initiations

    # Add these newly initiated molecules to the ones that were already initiated
    initiation_times.extend(final_new_initiations)
    start_label_pos_list.extend([1] * len(final_new_initiations))

    #Sample to always get 5000 molecules, this could have changed due to the initiation sampling
    desired_num = 5000
    # Ensure both lists are the same length before sampling
    assert len(initiation_times) == len(start_label_pos_list)
    # Sample unique indices if the combined pool is bigger than desired_num
    if len(initiation_times) > desired_num:
        sampled_indices = random.sample(range(len(initiation_times)), desired_num)
        initiation_times_sampled = [initiation_times[i] for i in sampled_indices]
        start_label_pos_list_sampled = [start_label_pos_list[i] for i in sampled_indices]
    else:
        initiation_times_sampled = initiation_times
        start_label_pos_list_sampled = start_label_pos_list

    initiation_times=initiation_times_sampled
    start_label_pos_list=start_label_pos_list_sampled


    #calculate stop position for each mRNA molecule
    for i in range(len(initiation_times)):
        initiation=initiation_times[i]
        initiation=initiation/60#Convert initiation time back to minutes
        start_label_pos=start_label_pos_list[i]

        # use function #4 to find the stop site nt coordinate
        stop_label_pos = findStopnt(df_nt, initiation, start_label_pos) 
        strand = df_nt.loc[0, 'strand']

        # construct the read sequence, using the coordinates identified from labeling, which was available for the entirety of transcription
        read = dna_sequence[0:int(stop_label_pos)]

        # convert positions
        converted_read, percent_sub, converted_positions, incorporated_positions = convert_base(
            read,
            nt_incorporation_rate,
            sub_rate_percent_range,
            int(start_label_pos),
            int(stop_label_pos),
            from_base,
            to_base
        )
        
        seq_err_read, percentage_seq_err, seq_err_positions = add_seq_errs(
            converted_read,
            list(map(float, args.seq_err.split(',')))
        )

        n = n + 1
        molecule_id = id_generator()

        #If an initiation time higher than the exp time was assigned (to avoid TSS pausing), then export correctly as the experiment time
        if initiation>args.experiment_time:
            initiation=args.experiment_time

        row_data = [
            initiation,
            molecule_id,
            strand,
            sub_rate_percent_range,
            start_label_pos,
            stop_label_pos,
            converted_positions,
            incorporated_positions,
            percentage_seq_err,
            seq_err_positions,
            seq_err_read
        ]
        data_list.append(row_data)
        # append only if at least one base was incorporated, if not, the molecule has no modified nucleotides and won't get enriched
        #if incorporated_positions:
        


    #Add background mRNAs that were made before the labeling (fully transcribed, unspliced, unlabeled)
    bkg_molecules=args.bkg_molecules
    if bkg_molecules > 0:
        bg_data_list = []
        output_filename_bg = str(args.o) + base_filename + "_background.tsv.gz"
    # Calculate how many nascent molecules to keep
        num_to_sample = round(len(data_list) * (1-bkg_molecules))
        data_list = random.sample(data_list, num_to_sample) #sample from the nascent molecules
        bkg_molecules_to_simulate=round(len(data_list)*bkg_molecules) #Get the number of background molecules
        full_sequence = ''.join(df['sequence'].astype(str)) #concatenate all the features for the transcript

        for _ in range(bkg_molecules_to_simulate):
            seq_err_read, percentage_seq_err, seq_err_positions = add_seq_errs(full_sequence, list(map(float, args.seq_err.split(',')))) #Add seq error
            #Add to the df
            molecule_id=id_generator()
            initiation=0
            strand=df_nt.loc[0, 'strand']
            sub_rate_percent_range=0
            start_label_pos=0
            stop_label_pos=0
            converted_positions='NA'
            incorporated_positions='NA'
            percentage_seq_err=percentage_seq_err
            seq_err_positions=seq_err_positions
            seq_err_read=seq_err_read
            
            bkg_row_data = [
                initiation, molecule_id, strand, sub_rate_percent_range, start_label_pos, stop_label_pos,
                converted_positions, incorporated_positions, percentage_seq_err, seq_err_positions, seq_err_read
            ]
            #Export background molecules to then incporate across simulated experiments
            bg_data_list.append(bkg_row_data)

        bg_data_list_for_export = pd.DataFrame(bg_data_list, columns=["initiation_time", 'molecule_id','strand',"sub_rate_percent_range", 'start_label_pos',
                                            'stop_label_pos','converted_positions','incorporated_positions', 'percentage_seq_err',
                                            'seq_err_positions', 'full_molecule_sequence'])
        # Save to a tab-separated file
        bg_data_list_for_export.to_csv(output_filename_bg, sep='\t', index=False, compression='gzip')


    # Convert to DataFrame with column names
    df_for_export = pd.DataFrame(data_list, columns=["initiation_time", 'molecule_id','strand',"sub_rate_percent_range", 'start_label_pos',
                                            'stop_label_pos','converted_positions','incorporated_positions', 'percentage_seq_err',
                                            'seq_err_positions', 'full_molecule_sequence'])


    # Save to a tab-separated file
    df_for_export.to_csv(output_filename, sep='\t', index=False, compression='gzip')


