# second set of functions - to be called after the first set

import pandas as pd
import random
import argparse
import math
import os
import string
import re

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
        
        # save the rate at that nucleotide for setting up the next start site
        #endsite_nt_rate = nt_traversal_time_df['rate_for_this_nt'].iloc[nt + 1] # for 1-based indexing
    
    #return(startsite_nt_coord, endsite_nt_coord, endsite_nt_rate)
    return(endsite_nt_coord)


# 5: Perform nucleotide conversions within the labeled region
#def convert_base(read, nt_incorporation_rate, nt_conversion_rate, start_site, stop_site, from_base, to_base):
    
    # this will read the sequence of the region created during the labeling period in question (5eU, 6sG, or 4sU) and put it into a substring
#    labeled_region = read[start_site:stop_site]
    
    # get all of the positions in the substring that correspond to the "from_base" character
#    pos = [pos for pos, char in enumerate(labeled_region) if char == from_base]
    
    # determine how many of the from_base bases will have modified nts incorporated 
#    num_bases_incorporated = int(math.ceil(len(pos) * nt_incorporation_rate))
    
    # collect all the positions that will have modified nucleotides incorporated
#    pos_conv = random.sample(pos, num_bases_incorporated)
    
    # convert string of the nt read to a list so nt characters can be changed
#    conv_read = list(read)
    
    # actually substitute the nts in the read using this loop
#    for i in pos_conv:
        
        # this needs to take into account the actual conversion rate for the given nucleotide (N/A for 5eU, 0.85-0.90 for 4sU, 0.35-0.55 for 6sG)
#        if random.random() < random.uniform(nt_conversion_rate[0], nt_conversion_rate[1]):
            
#            conv_read[start_site + i] = to_base
            
#    return ''.join(conv_read), nt_incorporation_rate, [start_site + i for i in pos_conv]

def convert_base(read, nt_incorporation_rate, nt_conversion_rate, start_site, stop_site, from_base, to_base):
    
    #Extract the labeled region of the read
    labeled_region = read[start_site:stop_site]
    
    # Find positions of all 'from_base' in the labeled region
    pos = [pos for pos, char in enumerate(labeled_region) if char == from_base]
    
    nt_incorporation_rate_to_use=random.uniform(nt_incorporation_rate[0], nt_incorporation_rate[1])
    #Determine how many will be incorporated
    num_bases_incorporated = int(math.ceil(len(pos) * nt_incorporation_rate_to_use))
    
    #Randomly choose which ones are incorporated
    pos_conv = random.sample(pos, num_bases_incorporated)
    
    #Determine what percent of those will actually be converted
    conversion_fraction = random.uniform(nt_conversion_rate[0], nt_conversion_rate[1])
    num_bases_converted = int(math.floor(len(pos_conv) * conversion_fraction))
    
    #Choose randomly which of the incorporated bases are converted
    pos_to_convert = random.sample(pos_conv, num_bases_converted)
    
    #Convert to list and apply the conversion
    conv_read = list(read)
    for i in pos_to_convert:
        conv_read[start_site + i] = to_base #->>>>>WHY IS THIS ADDING TO STARTSITE?
    
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


# 7: create the output name for the gtf to feed into the chopper script
def createOutputName(gtf_file, label_time, bkgd_perc_range, u2c_perc_range, num_reads):
    
    gtf_name = gtf_file.strip(".gtf")
    
    output_file_name = str(gtf_name) + "_l" + str(label_time) + "_b" + str(bkgd_perc_range) + "_p" + str(u2c_perc_range) + "_reads" + str(num_reads) + ".csv.gz"

#8: add molecule id (this will also be the read name for long reads)
def id_generator(size=12, chars=string.ascii_uppercase + string.digits):
	return "@" + ''.join(random.choice(chars) for _ in range(size - 1))   

#def str_to_float_list(s):
#    return [float(x) for x in s.split(',')]

# now we call the main script
if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument('--region_file', default='tmp.gtf', help='Input region specific elongation rate GTF file name')
    parser.add_argument('--nt_file', default='tmp.gtf', help='Input ntTraversalTime GTF file name')
    parser.add_argument('--l', type=int, default=15, help='Labeling time in minutes')
    parser.add_argument('--seq_err', default='0.0001,0.0002', help='Illumina sequencing error rate range')
    parser.add_argument('--nt_inc_rate',type=str,default='0.09,0.1',help='Comma-separated range of nucleotide incorporation proportions')
    parser.add_argument('--subs_rate',type=str,default='0.95,1',help='Comma-separated range of nucleotide substitution proportions')
    parser.add_argument('--sub_type', type=str, default='T,T', help='Specify the expected conversion for the nucleotide. For example T,C for 4sU. Default is no substitution.', required=False)
    #parser.add_argument('--seq_depth', type=int, default=20000000, help='Total library sequencing depth')
    parser.add_argument('--bkg_molecules', type=float,default=0,help='Proportion of molecules that are background (0 or 1)')
    parser.add_argument('--path_to_tpm', type=str,help='Full path to the TPM per gene file')
    parser.add_argument('--o', type=str, default='./', help='output path')
    parser.add_argument('--treatment', type=str, default='no DRB', help='usage of DRB or not')
    args = parser.parse_args()


    # set variables for arguments needed for functions 4-7
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
    output_filename = str(args.o) + base_filename + "_monolabeled.tsv.gz"

    # open the output file for editing
    output_file = open(output_filename, 'w')

    # initialize the read number counter
    n=0

    # initialize the master data list that will be output at the end
    data_list = []

    #Get expected subs
    from_base, to_base = args.sub_type.split(',')


    max_transcript_length = df["region_end_coord"].max()

    # Estimate molecule count for this single transcript based on TPM
    #Calculate an approximate max transcribed length using the mean elong rate across nucleotides
    
    #Simulate 2000 mRNA molecules. This is based on the mean number of molecules needed to get a TPM of 200 (upper limit)
    molecule_count = 2000


    #Generate initiation times for generating reads, one for each read spanning the entire labeling experiment (labeling time x 3)
    # initiation time is treated as the time remaining in the experiment, for example if initiation = 2, then there is 2 minutes remaining in the experiment
    initiation_times=[]
    for i in range(0,molecule_count):
        initiation_temp=random.uniform(0,args.l) 
        initiation_times.append(initiation_temp)
    

    # the underscore here means the value has no import on the function (here it is the read number)
    for _ in range(molecule_count):
        # grab the nth initiation time
        test = initiation_times[n]
        initiation = test

        #If DRB treated labeling starts in pos 1 of the transcript, if not it's random
        #start_label_pos = 1 if args.treatment == "DRB" else random.randint(1, max_transcript_length)
        start_label_pos = 1 if args.treatment == "DRB" else random.randint(-10000, max_transcript_length) #-10000 is an imaginary position for Pol IIs that will start transcription after the start of the experiment

        #initiation = args.l if args.treatment == "no DRB" else initiation == initiation
        initiation = test if args.treatment != "no DRB" else args.l


        #Add position for no DRB molecules that initiate after the experiment
        if start_label_pos<1:#If pol II hasn't started yet
            time_needed_for_tss=abs(start_label_pos)/2000 #Time required to get to the TSS
            initiation=initiation-time_needed_for_tss #update time left 
            start_label_pos=1#Update the position to the TSS
            if initiation<=0:
                initiation=0#If Pol II never gets to the TSS then it doesn't have time to elongate


        
        # use function #4 to find the stop site nt coordinate for 5eU labeling (the first 5 minutes)
        # use initiation as the labeling duration, as 5eU will be available for the entirety of the experiment (if it doesn't reach the end of the gene, then that transcript wouldn't be fully transcribed)
        stop_label_pos = findStopnt(df_nt, initiation, start_label_pos) 
        strand=df_nt.loc[0, 'strand']

        # construct the read sequence, using the coordinates identified from 5eU labeling, which was available for the entirety of transcription
        read = dna_sequence[1:int(stop_label_pos)]

        # convert positions
        converted_read, percent_sub, converted_positions = convert_base(read, nt_incorporation_rate, sub_rate_percent_range, int(start_label_pos), int(stop_label_pos),from_base,to_base)
        
        #bg_mutated_read, percentage_bg_mutations, mutated_bg_positions = mutate_bg(converted_read, list(map(float, args.b.split(','))), int(starts_5eU_pos), int(stop_site_5eU))
        seq_err_read, percentage_seq_err, seq_err_positions = add_seq_errs(converted_read, list(map(float, args.seq_err.split(','))))

        n = n + 1

        molecule_id=id_generator()


        row_data = [initiation, molecule_id,strand,sub_rate_percent_range, start_label_pos, stop_label_pos, converted_positions, percentage_seq_err,seq_err_positions,seq_err_read]
        #row_data.columns = ["initiation_time", "sub_rate_percent_range", 'start_label_pos','stop_label_pos', 'converted_positions', 'percentage_seq_err','seq_err_positions','full_molecule_sequence']

        # Add the row_data to the data_list

        data_list.append(row_data)

    
    #Add background mRNAs that were made before the labeling (fully transcribed, unspliced, unlabeled)
    bkg_molecules=args.bkg_molecules
    if bkg_molecules > 0:
    # Calculate how many nascent molecules to keep
        num_to_sample = round(len(row_data) * (1-bkg_molecules))
        row_data = random.sample(row_data, num_to_sample) #sample from the nascent molecules
        bkg_molecules_to_simulate=round(len(row_data)*bkg_molecules) #Get the number of background molecules
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
            percentage_seq_err=percentage_seq_err
            seq_err_positions=seq_err_positions
            seq_err_read=seq_err_read
            data_list.append(row_data)

    # Convert to DataFrame with column names
    df_for_export = pd.DataFrame(data_list, columns=["initiation_time", 'molecule_id','strand',"sub_rate_percent_range", 'start_label_pos',
                                      'stop_label_pos','converted_positions', 'percentage_seq_err',
                                      'seq_err_positions', 'full_molecule_sequence'])


    # Save to a tab-separated file
    df_for_export.to_csv(output_filename, sep='\t', index=False, compression='gzip')







