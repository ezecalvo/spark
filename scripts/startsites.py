import sys
import pdb
import os
import subprocess
import argparse
import collections
from os.path import isfile, join
import gzip
import getopt
import re
import pysam

def readBam(bamfile):
    bamhere = pysam.Samfile(bamfile, 'rb')
    return(bamhere)

def startsites_single(bamhere, bamfile, readtype, includesplit=False):
    bamname = bamfile[:-4]
    if includesplit == False:
        outfile = gzip.open(bamname +'_startsites_read1.bed.gz', 'wt')
        for read in bamhere:
            if 'N' not in read.cigarstring:
                start = int(read.pos)
                strand = '+'
                if read.is_reverse:
                    start = int(read.aend)
                    strand = '-'
                end = int(start) + 1
                outfile.write(str(bamhere.getrname(read.rname)) +'\t'+ str(start) + '\t'+ str(end) + '\t'+ str(read.qname) +'\t.\t'+ str(strand) + '\n')
    if includesplit == True:
        outfile = gzip.open(bamname +'_startsites_read1_withJunctions.bed.gz', 'wt')
        for read in bamhere:
            start = int(read.pos)
            strand = '+'
            if read.is_reverse:
                start = int(read.aend)
                strand = '-'
            end = int(start) + 1
            outfile.write(str(bamhere.getrname(read.rname)) +'\t'+ str(start) + '\t'+ str(end) + '\t'+ str(read.qname) +'\t.\t'+ str(strand) + '\n')
    bamhere.close()
    outfile.close()

def startsites_paired(bamhere, bamfile, readtype, includesplit=False):
    bamname = bamfile[:-4]
    if includesplit == False:
        outfile1 = gzip.open(bamname +'_startsites_read1.bed.gz', 'wt')
        outfile2 = gzip.open(bamname +'_startsites_read2.bed.gz', 'wt')
        for read in bamhere:
            if read.is_paired and read.is_proper_pair and 'N' not in read.cigarstring:
                start = int(read.pos)
                strand = "+"
                if read.is_reverse:
                    start = int(read.aend)
                    strand = '-'
                end = int(start) + 1
                if read.is_read1:
                    outfile1.write(str(bamhere.getrname(read.rname)) +'\t'+ str(start) + '\t'+ str(end) + '\t'+ str(read.qname) +'\t.\t'+ str(strand) + '\n')
                if not read.is_read1:
                    outfile2.write(str(bamhere.getrname(read.rname)) +'\t'+ str(start) + '\t'+ str(end) + '\t'+ str(read.qname) +'\t.\t'+ str(strand) + '\n')
    if includesplit == True:
        outfile1 = gzip.open(bamname +'_startsites_withJunctions_read1.bed.gz', 'wt')
        outfile2 = gzip.open(bamname +'_startsites_withJunctions_read2.bed.gz', 'wt')
        for read in bamhere:
            if read.is_paired and read.is_proper_pair:
                start = int(read.pos)
                strand = "+"
                if read.is_reverse:
                    start = int(read.aend)
                    strand = '-'
                end = int(start) + 1
                if read.is_read1:
                    outfile1.write(str(bamhere.getrname(read.rname)) +'\t'+ str(start) + '\t'+ str(end) + '\t'+ str(read.qname) +'\t.\t'+ str(strand) + '\n')
                if not read.is_read1:
                    outfile2.write(str(bamhere.getrname(read.rname)) +'\t'+ str(start) + '\t'+ str(end) + '\t'+ str(read.qname) +'\t.\t'+ str(strand) + '\n')
    bamhere.close()
    outfile1.close()
    outfile2.close()

def startsites_pairedFirst(bamhere, bamfile, readtype, includesplit=False):
    bamname = bamfile[:-4]
    if includesplit == False:
        outfile1 = gzip.open(bamname +'_startsites_firstInPair_read1.bed.gz', 'wt')
        outfile2 = gzip.open(bamname +'_startsites_firstInPair_read2.bed.gz', 'wt')
        for read in bamhere:
            if read.is_paired and read.is_proper_pair and read.is_read1 and 'N' not in read.cigarstring:
                start = int(read.pos)
                strand = "+"
                if read.is_reverse:
                    start = int(read.aend)
                    strand = '-'
                end = int(start) + 1
                if read.is_read1:
                    outfile1.write(str(bamhere.getrname(read.rname)) +'\t'+ str(start) + '\t'+ str(end) + '\t'+ str(read.qname) +'\t.\t'+ str(strand) + '\n')
                if not read.is_read1:
                    outfile2.write(str(bamhere.getrname(read.rname)) +'\t'+ str(start) + '\t'+ str(end) + '\t'+ str(read.qname) +'\t.\t'+ str(strand) + '\n')                
    if includesplit == True:
        outfile1 = gzip.open(bamname +'_startsites_firstInPair_withJunctions_read1.bed.gz', 'wt')
        outfile2 = gzip.open(bamname +'_startsites_firstInPair_withJunctions_read2.bed.gz', 'wt')
        for read in bamhere:
            if read.is_paired and read.is_proper_pair and read.is_read1:
                start = int(read.pos)
                strand = "+"
                if read.is_reverse:
                    start = int(read.aend)
                    strand = '-'
                end = int(start) + 1
                if read.is_read1:
                    outfile1.write(str(bamhere.getrname(read.rname)) +'\t'+ str(start) + '\t'+ str(end) + '\t'+ str(read.qname) +'\t.\t'+ str(strand) + '\n')
                if not read.is_read1:
                    outfile2.write(str(bamhere.getrname(read.rname)) +'\t'+ str(start) + '\t'+ str(end) + '\t'+ str(read.qname) +'\t.\t'+ str(strand) + '\n')
    bamhere.close()
    outfile1.close()
    outfile2.close()

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--bamfile',type=str,help='bam file', required=True)
    #parser.add_argument('--outbed',type=str,help='bed file for output', required = True)
    parser.add_argument('--readtype', type=str, default = 'paired', choices=['single','paired'], help='single or paired-end reads (default = paired)', required = True)
    parser.add_argument('--includesplit', action='store_true', default=False, help='include split reads (default = F)', required=False)
    parser.add_argument('--onlyfirstread',action='store_true', default=False, help='only first read in pair (default = F)', required=False)
    args = parser.parse_args()

    ## read in bam files
    bamhere = readBam(args.bamfile)

    if args.readtype == 'single' and args.split == False:
        startsites_single(bamhere, args.bamfile, args.readtype, includesplit=False)
    if args.readtype == 'single' and args.split == True:
        startsites_single(bamhere, args.bamfile, args.readtype, includesplit=True)

    if args.readtype == 'paired' and args.onlyfirstread == False and args.split == False:
        startsites_paired(bamhere, args.bamfile, args.readtype, includesplit=False)
    if args.readtype == 'paired' and args.onlyfirstread == False and args.split == True:
        startsites_paired(bamhere, args.bamfile, args.readtype, includesplit=True)

    if args.readtype == 'paired' and args.onlyfirstread == True and args.split == False:
        startsites_pairedFirst(bamhere, args.bamfile, args.readtype, includesplit=False)
    if args.readtype == 'paired' and args.onlyfirstread == True and args.split == True:
        startsites_pairedFirst(bamhere, args.bamfile, args.readtype, includesplit=True)
