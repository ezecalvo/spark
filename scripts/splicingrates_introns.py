# script to get intron information out of a gtf file
# input is a gtf file of annotations
# output is one bed6 file per desired intron type


#ECR 3-9-26: Made some changes to accurately get all and alternative exons without errors. getAIs was the same logic as getCIs



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
    if gtf.endswith(".gz"): 
        gtffile = gzip.open(gtf, 'rt') 
    else: 
        gtffile = open(gtf, 'r')
        
    genedict = collections.defaultdict(lambda: collections.defaultdict(dict))
    gene = ""
    transcript = ""
    
    for line in gtffile:
        if line.startswith("#"): continue
        gtfline = line.strip().split('\t')
        
        # Skip incomplete lines
        if len(gtfline) < 9: continue
        
        feature = gtfline[2]
        if feature not in ["gene", "transcript", "exon", "start_codon", "stop_codon"]:
            continue
            
        # ECR 3-9-26: the previous version was not spliting the 9th column correctly and transcripts were not getting merge properly 
        attr_dict = {}
        for attr in gtfline[8].split(";"):
            attr = attr.strip()
            if not attr: continue
            parts = attr.split(" ")
            if len(parts) >= 2:
                attr_dict[parts[0]] = parts[1].strip('"')

        if feature == "gene":
            gene_id = attr_dict.get("gene_id", "unknown_gene")
            gene = f"{gene_id}:{gtfline[0]}:{gtfline[3]}-{gtfline[4]}:{gtfline[6]}"
            
        elif feature == "transcript":
            transcript_id = attr_dict.get("transcript_id", "unknown_transcript")
            transcript = f"{transcript_id}:{gtfline[0]}:{gtfline[3]}-{gtfline[4]}"
            genedict[gene][transcript]['exon_starts'] = []
            genedict[gene][transcript]['exon_ends'] = []
            
        elif feature == "exon":
            genedict[gene][transcript]['exon_starts'].append(gtfline[3])
            genedict[gene][transcript]['exon_ends'].append(gtfline[4])
            
        elif feature == "start_codon":
            genedict[gene][transcript]['start'] = [gtfline[3], gtfline[4]]
            
        elif feature == "stop_codon":
            genedict[gene][transcript]['end'] = [gtfline[3], gtfline[4]]
            
    print(str(len(genedict)) + " genes read.")
    gtffile.close()
    return genedict



# all introns
def getIntrons(genedict,outname):
    print("Reading introns.")
    outfile=open(outname +".allintrons.bed",'w')
    outfiledist = open(outname +".allintrons.threedist", 'w')
    outfiledist.write("name\tdistance\n")
    intronlen=0
    for gene in genedict:
        introns = set() 
        strand=gene.split(":")[3]
        threeend = int(gene.split(":")[2].split("-")[1])
        if strand == "-": threeend = int(gene.split(":")[2].split("-")[0])
        
        for transcript in genedict[gene]:
            # ECR 3-9-26: Convert to integers before sorting
            exonstarts = sorted([int(x) for x in genedict[gene][transcript]['exon_starts']])
            exonends = sorted([int(x) for x in genedict[gene][transcript]['exon_ends']])
            for i in range(0,(len(exonends)-1)):
                introns.add(str(exonends[i]) +"-"+ str(exonstarts[i+1]))
                
        intronlen += len(introns)
        for intron in introns:
            name = str(gene.split(":")[0]) +':'+ str(gene.split(":")[1]) +':'+ str(intron) +':'+ str(gene.split(":")[3])
            outfile.write(str(gene.split(":")[1]) +"\t"+ str(intron.split("-")[0]) +"\t"+ str(intron.split("-")[1]) +"\t"+ str(name) +"\tintron\t"+ str(gene.split(":")[3]) +"\n")#ECR 3-9-26: Changed 'one' to 'intron'
            
            # write file with distance to 3' gene end (for kinetics scripts)
            threedist = threeend - int(intron.split("-")[1])
            if strand == "-": threedist = int(intron.split("-")[0]) - threeend
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
        consts = None 
        strand=gene.split(":")[3]
        threeend = int(gene.split(":")[2].split("-")[1])
        if strand == "-": threeend = int(gene.split(":")[2].split("-")[0])
        
        for transcript in genedict[gene]:
            intronnames = set()
            exonstarts = sorted([int(x) for x in genedict[gene][transcript]['exon_starts']])
            exonends = sorted([int(x) for x in genedict[gene][transcript]['exon_ends']])
            for i in range(0,(len(exonends)-1)): 
                intronnames.add(str(exonends[i]) +"-"+ str(exonstarts[i+1]))
            
            #ECR 3-9-26: Use set intersections to find constitutive introns safely
            if consts is None:
                consts = intronnames
            else:
                consts = consts.intersection(intronnames)
                
            if len(consts) == 0: break # #ECR 3-9-26: if no constitutive introns left, break
            
        if consts is None: consts = set() 
        constlenall += len(consts)
        
        for intron in consts:
            name = str(gene.split(":")[0]) +':'+ str(gene.split(":")[1]) +':'+ str(intron) +':'+ str(gene.split(":")[3])
            outfile.write(str(gene.split(":")[1]) +"\t"+ str(intron.split("-")[0]) +"\t"+ str(intron.split("-")[1]) +"\t"+ str(name) +"\tintron\t"+ str(gene.split(":")[3]) +"\n")
            
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
    altlen=0
    for gene in genedict:
        all_introns = set()
        consts = None
        strand=gene.split(":")[3]
        threeend = int(gene.split(":")[2].split("-")[1])
        if strand == "-": threeend = int(gene.split(":")[2].split("-")[0])
        
        for transcript in genedict[gene]:
            intronnames = set()
            exonstarts = sorted([int(x) for x in genedict[gene][transcript]['exon_starts']])
            exonends = sorted([int(x) for x in genedict[gene][transcript]['exon_ends']])
            for i in range(0,(len(exonends)-1)): 
                intronnames.add(str(exonends[i]) +"-"+ str(exonstarts[i+1]))
            
            all_introns.update(intronnames)
            
            if consts is None:
                consts = set(intronnames)
            else:
                consts = consts.intersection(intronnames)
        
        if consts is None: consts = set()
        
        # ECR 3-9-26: alternative introns where not being parsed correctly. Now alternative is all exons except constitutive
        alt_introns = all_introns - consts
        altlen += len(alt_introns)
        
        for intron in alt_introns:
            name = str(gene.split(":")[0]) +':'+ str(gene.split(":")[1]) +':'+ str(intron) +':'+ str(gene.split(":")[3])
            outfile.write(str(gene.split(":")[1]) +"\t"+ str(intron.split("-")[0]) +"\t"+ str(intron.split("-")[1]) +"\t"+ str(name) +"\tintron\t"+ str(gene.split(":")[3]) +"\n")
            
            threedist = threeend - int(intron.split("-")[1])
            if strand == "-": threedist = int(intron.split("-")[0]) - threeend
            outfiledist.write(str(name) +"\t"+ str(threedist) +"\n")
            
    print("... " +str(altlen)+ " alternative introns found.")
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
        
        
        
        











