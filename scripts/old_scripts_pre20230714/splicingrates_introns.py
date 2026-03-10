# script to get intron information out of a gtf file
# input is a gtf file of annotations
# output is one bed6 file per desired intron type

# Usage: python splicingrates_introns.py --help

import sys
import pdb
import os
import subprocess
import argparse
import collections
from os.path import isfile, join
import gzip
import getopt

# read in gtf to create dictionary with one entry per gene
def getGenes(gtf):
    if gtf[-2:] == "gz": gtffile=gzip.open(gtf)
    else: gtffile=open(gtf)
    genedict=collections.defaultdict(lambda:collections.defaultdict(dict))
    while('TRUE'):
        gtfline=gtffile.readline().split()
        if not gtfline: break
        if gtfline[0][:1] == "#": continue
        if gtfline[2] == "gene":
            gene=gtfline[9].split("\"")[1]+":"+gtfline[0]+":"+gtfline[3]+"-"+gtfline[4]+":"+gtfline[6]
        if gtfline[2] == "transcript":
            transcript=gtfline[11].split("\"")[1]+":"+gtfline[0]+":"+gtfline[3]+"-"+gtfline[4]
            genedict[gene][transcript]['exon_starts'] = []
            genedict[gene][transcript]['exon_ends'] = []
        if gtfline[2] == "exon":
            genedict[gene][transcript]['exon_starts'].append(gtfline[3])
            genedict[gene][transcript]['exon_ends'].append(gtfline[4])
        if gtfline[2] == "start_codon":
            genedict[gene][transcript]['start'] = [gtfline[3],gtfline[4]]
        if gtfline[2] == "stop_codon":
            genedict[gene][transcript]['end'] = [gtfline[3],gtfline[4]]
    print(str(len(genedict)) +" genes read.")
    return(genedict)

# all introns
def getIntrons(genedict,outname):
    print("Reading introns.")
    outfile=open(outname +".allintrons.bed",'w')
    outfiledist = open(outname +".allintrons.threedist", 'w')
    outfiledist.write("name\tdistance\n")
    intronlen=0
    for gene in genedict:
        introns=[]
        strand=gene.split(":")[3]
        threeend = int(gene.split(":")[2].split("-")[1])
        if strand == "-": threeend = int(gene.split(":")[2].split("-")[0])
        for transcript in genedict[gene]:
            exonstarts = sorted(genedict[gene][transcript]['exon_starts'])
            exonends = sorted(genedict[gene][transcript]['exon_ends'])
            for i in range(0,(len(exonends)-1)):
                introns.append(exonends[i] +"-"+ exonstarts[i+1])
        introns2=list(set(introns))
        intronlen += len(introns2)
        for intron in introns2:
            name = str(gene.split(":")[0]) +':'+ str(gene.split(":")[1]) +':'+ str(intron) +':'+ str(gene.split(":")[3])
            outfile.write(str(gene.split(":")[1]) +"\t"+ str(one.split("-")[0]) +"\t"+ str(one.split("-")[1]) +"\t"+ str(name) +"\tintron\t"+ str(gene.split(":")[3]) +"\n")
            # write file with distance to 3' gene end (for kinetics scripts)
            threedist = threeend - int(one.split("-")[1])
            if strand == "-": threedist = int(one.split("-")[0]) - threeend
            outfiledist.write(str(name) +"\t"+ str(threedist) +"\n")
    print("... " +str(intronlen)+ " introns found.")
    outfile.close()
    outfiledist.close()

# constitutive introns
def getCIs(genedict,outname):
    print("Reading constitutive introns.")
    outfile=open(outname +".constitutive_introns.bed",'w')
    outfiledist = open(outname +".constitutive_introns.threedist", 'w')
    outfiledist.write("name\tdistance\n")
    constlenall=0
    for gene in genedict:
        consts=[]
        strand=gene.split(":")[3]
        threeend = int(gene.split(":")[2].split("-")[1])
        if strand == "-": threeend = int(gene.split(":")[2].split("-")[0])
        for transcript in genedict[gene]:
            intronnames = []
            exonstarts = sorted(genedict[gene][transcript]['exon_starts'])
            exonends = sorted(genedict[gene][transcript]['exon_ends'])
            for i in range(0,(len(exonends)-1)): 
                intronnames.append(exonends[i] +"-"+ exonstarts[i+1])
            if not consts: consts = intronnames
            for intron in consts:
            	if intron not in intronnames: consts.remove(intron)
            constlen = len(consts)
            if constlen == 0: break
        constlenall += len(consts)
        for intron in consts:
            name = str(gene.split(":")[0]) +':'+ str(gene.split(":")[1]) +':'+ str(intron) +':'+ str(gene.split(":")[3])
            outfile.write(str(gene.split(":")[1]) +"\t"+ str(intron.split("-")[0]) +"\t"+ str(intron.split("-")[1]) +"\t"+ str(name) +"\tintron\t"+ str(gene.split(":")[3]) +"\n")
            # write file with distance to 3' gene end (for kinetics scripts)
            threedist = threeend - int(intron.split("-")[1])
            if strand == "-": threedist = int(intron.split("-")[0]) - threeend
            outfiledist.write(str(name) +"\t"+ str(threedist) +"\n")
    print("... " +str(constlenall)+ " constitutive introns found.")
    outfile.close()
    outfiledist.close()

# alternative introns
def getAIs(genedict,outname):
    print("Reading alternative introns.")
    outfile=open(outname +".alternative_introns.bed",'w')
    outfiledist = open(outname +".alternative_introns.threedist", 'w')
    outfiledist.write("name\tdistance\n")
    constlen=0
    for gene in genedict:
        consts=[]
        strand=gene.split(":")[3]
        threeend = int(gene.split(":")[2].split("-")[1])
        if strand == "-": threeend = int(gene.split(":")[2].split("-")[0])
        for transcript in genedict[gene]:
            intronnames = []
            exonstarts = sorted(genedict[gene][transcript]['exon_starts'])
            exonends = sorted(genedict[gene][transcript]['exon_ends'])
            for i in range(0,(len(exonends)-1)): 
                intronnames.append(exonends[i] +"-"+ exonstarts[i+1])
            if not consts: consts = intronnames
            for intron in consts:
            	if intron not in intronnames: consts.remove(intron)
        constlen += len(consts)
        for intron in consts:
            name = str(gene.split(":")[0]) +':'+ str(gene.split(":")[1]) +':'+ str(intron) +':'+ str(gene.split(":")[3])
            outfile.write(str(gene.split(":")[1]) +"\t"+ str(intron.split("-")[0]) +"\t"+ str(intron.split("-")[1]) +"\t"+ str(name) +"\tintron\t"+ str(gene.split(":")[3]) +"\n")
            # write file with distance to 3' gene end (for kinetics scripts)
            threedist = threeend - int(intron.split("-")[1])
            if strand == "-": threedist = int(intron.split("-")[0]) - threeend
            outfiledist.write(str(name) +"\t"+ str(threedist) +"\n")
    print("... " +str(constlen)+ " constitutive introns found.")
    outfile.close()
    outfiledist.close()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                     description='Identify introns for estimating splicing rates.', add_help = True)
    parser.add_argument('--gtf', type=str, metavar='X.gtf', help='gtf file from which to get annotations (full path)', required=True)
    parser.add_argument('--introns',type=str, default = 'constitutive', choices=['allintrons','alternative','constitutive'], help='comma separated list of types of annotation to output into a file')
    parser.add_argument('--outname',type=str, metavar='', help='basename of output file (full path)', required=True)
    args = parser.parse_args()

    genedict=getGenes(args.gtf)

    #pdb.set_trace()    

    regions=args.introns
    # so far, types include: 5/3UTRs, 5/3splicesites, coding exons, constitutive exons, introns
    for type in regions.split(","):
#        print str(type)
        if type == "allintrons": getIntrons(genedict,args.outname)
        if type == "constitutive": getCIs(genedict, args.outname)
        if type == "alternative": getAIs(genedict, args.outname)
