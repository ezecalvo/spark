# script to calculate splicing rates for designated introns
# Usage: python splicingrates_model.py --help

import sys
import pdb
import os
import subprocess
import argparse
import collections
from os.path import isfile, join
import gzip
import getopt
import math
import pandas as pd

def readsPerReplicate(basename, times, rep):
	suffix = 'm_rep' + str(rep) 
	name_t1 = basename +'_'+ times[0]+ suffix + '_junctionCombo.coverage.gz'
	df_rep = pd.read_csv(name_t1, sep='\t', compression="gzip").set_index('intron')
	#df_rep.columns = ['intron', times[0] + suffix +'_ie', times[0] + suffix + '_ee']
	df_rep.columns = [times[0] +'m_ie', times[0] + 'm_ee']
	df_rep[times[0]+'m_ratio'] = df_rep[times[0]+'m_ie'] / df_rep[times[0]+'m_ee']
	if len(times) > 1:
		for t in range(1, len(times)):
			name_t = basename +'_'+ times[t]+ suffix + '_junctionCombo.coverage.gz'
			df_here = pd.read_csv(name_t, sep='\t', compression="gzip").set_index('intron')
			df_rep[times[t]+'m_ie'] = df_here['ie_count']
			df_rep[times[t]+'m_ee'] = df_here['ee_count']
			df_rep[times[t]+'m_ratio'] = df_here[times[t]+'m_ie'] / df_here[times[t]+'m_ee']
	return(df_rep)

def SummedReads(basename, times, reps):	
	## go through timepoint 1
	name_t1_r1 = times[0] + 'm_rep1'
	file_t1_r1 = basename +'_'+ name_t1_r1 + '_junctionCombo.coverage.gz'
	df_time = pd.read_csv(file_t1_r1, sep='\t', compression="gzip").set_index('intron')
	#df_time.columns = ['intron', name_t1_r1+'_ie', name_t1_r1+'_ee']
	df_time.columns = [name_t1_r1+'_ie', name_t1_r1+'_ee']
	if reps > 1:
		for r in range(1, reps):
			rep = times[0] + 'm_rep' + str(r+1)
			name_r = basename +'_'+ rep + '_junctionCombo.coverage.gz'
			df_here = pd.read_csv(name_r, sep='\t', compression="gzip").set_index('intron')
			df_time[rep+'_ie'] = df_here['ie_count']
			df_time[rep+'_ee'] = df_here['ee_count']
	df_time[times[0]+'m_ie'] = df_time.filter(like='_ie').sum()
	df_time[times[0]+'m_ee'] = df_time.filter(like='_ee').sum()
	df_time[times[0]+'m_ratio'] = df_time[times[0]+'m_ie'] / df_time[times[0]+'m_ee']
	# initialize new dataframe with time 1
	df_sum = pd.DataFrame(df_time, columns=["intron", times[0]+'m_ie', times[0]+'m_ee', times[0]+'m_ratio'])
	if len(times) > 1:
		for t in range(1, len(times)):
			name_rep1 = basename +'_'+ times[t] +'_m_rep1_junctionCombo.coverage.gz'
			df_time = pd.read_csv(name_rep1, sep='\t', compression="gzip").set_index('intron')
			#df_time.columns = ['intron', times[t]+'m_rep1_ie', times[t]+'m_rep1_ee']
			df_time.columns = [times[t]+'m_rep1_ie', times[t]+'m_rep1_ee']
			if reps > 1:
				for r in range(1, reps):
					rep = times[t] +'_m_rep'+ str(r+1)
					name_r = basename +'_'+ rep +'_junctionCombo.coverage.gz'
					df_here = pd.read_csv(name_r, sep='\t', compression="gzip").set_index('intron')
					df_time[rep+'_ie'] = df_here['ie_count']
					df_time[rep+'_ee'] = df_here['ee_count']
			df_sum[times[t]+'m_ie'] = df_time.filter(like='_ie').sum()
			df_sum[times[t]+'m_ee'] = df_time.filter(like='_ee').sum()
			df_sum[times[t]+'m_ratio'] = df_sum[times[t]+'m_ie'] / df_sum[times[t]+'m_ee']
	return(df_sum)

# read in 3' dist
def getThreeDist(threedistfile):
	df_threedist = pd.read_csv(threedistfile, sep='\t').set_index('name')
	return(df_threedist)


def calculateDprime(df_threelength, times, txnrate):
	# total.juncratio.data$D_prime <- (total.juncratio.data$threelength) + (minute * txnrate)
	for t in times:
		min = int(t)
		df_threelength[t+'m_Dprime'] = df_threelength['distance'] + (min * txnrate)
	return(df_threelength)

def calculateRprime(df_juncs, times, txnrate):
	#total.juncratio.data$R_prime <- (total.juncratio.data$D_prime*log(2))/(txnrate*((1/total.juncratio.data$ratio) + 1))
	for t in times:
		df_juncs[t+'m_Rprime'] = (df_juncs[t+'m_Dprime'] * math.log(2, 10)) / (txnrate * ((1 / df_juncs[t+'m_ratio']) + 1))
	return(df_juncs)

def estimateSplicingRates(halflifefile, txnrate):
	#print('CALCULATING SPLICING RATES')
	command = "Rscript splicingRatesModel.R " + halflifefile +" "+ str(txnrate)
	subprocess.call(command, shell=True)


if __name__ == '__main__':
	parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                     description='Estimate splicing rates using ie and ee junction reads.', add_help = True)
	### read counts
	group_reads = parser.add_argument_group('read information')
	group_reads.add_argument('--basename', type=str, metavar='x_Xm_repX', help="basename of combo files, including path. expected to have '_Xm_repX' prefix", required=True)
	group_reads.add_argument('--timepoints', type=str, metavar='t1,t2,...tn', help="comma separated list of 4sU labeling timepoints (minutes)", required=True)
	group_reads.add_argument('--replicates', type=int, default=1, metavar='', help="number of replicates for each timepoint")
	group_reads.add_argument('--summed', action='store_true', help='Estimate splicing rates with summed junction reads (across replicates)', required=False)
	### model assumptions
	group_model = parser.add_argument_group('model assumptions')
	group_model.add_argument('--outname',type=str, help='name for final splicing rate output file (full path)', required=True)
	group_model.add_argument('--txnrate', type=int, metavar='', default = 1500, help="transcription elongation rate in nt/min (default = 1500 nt/min)", required=False)
	group_model.add_argument('--threedist', type=str, metavar='x_threedist.txt', help="file with 3' distances for each intron from splicinrates_introns.py (full path)", required=True)

	args = parser.parse_args()

	print("Modeling Splicing Rates...")
	# parse timepoints
	times = args.timepoints.split(',')

	# read in 3' dist
	print("... getting 3' distances")
	df_threedist = getThreeDist(args.threedist)
	# calculate Dprime
	print("... calculating D'")
	df_threedist = calculateDprime(df_threedist, times, args.txnrate)

	# calculate rates for each replicates (across timepoints)
	print("... calculating rates per replicates")
	for r in range(0,args.replicates):
		rep = r+1
		print("...... rep{0}".format(rep))
		print("......... getting junction counts")
		df_rep = readsPerReplicate(args.basename, times, rep)
		# add Dprime values
		print("......... modeling half-lives")    	
		df_rep_Dprime = pd.concat([df_rep,df_threedist], axis=1)
		# calculate R'
		df_rep_Dprime_Rprime = calculateRprime(df_rep_Dprime, times, args.txnrate)
		# write outfile with junctions, D', and R'
		halflifefile = args.outname + '_rep' + str(r) + '.halflives' 
		df_rep_Dprime_Rprime.to_csv(halflifefile, sep='\t', index=True)
		estimateSplicingRates(halflifefile, args.txnrate)

	#pdb.set_trace()

	# calculate rates for summed replicate reads across timepoints
	if args.summed:
		print("... calculating rates for summed counts across replicates")
		print("...... getting junction counts")
		df_sum = SummedReads(args.basename, times, args.replicates)
		# add Dprime values
		print("......... modeling half-lives")
		df_sum_Dprime = pd.concat([df_sum,df_threedist], axis=1)
		# calculate R'
		df_sum_Dprime_Rprime = calculateRprime(df_sum_Dprime, times, args.txnrate)
		# write outfile with junctions, D', and R'
		halflifefile = args.outname + '_summed.halflives' 
		df_sum_Dprime_Rprime.to_csv(halflifefile, sep='\t', index=True)
		estimateSplicingRates(halflifefile, args.txnrate)
