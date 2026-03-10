# script to get particular annotation information out of a gtf file
# input is a gtf file of annotations
# output is one bed6 file per desired annotation type

# Usage: python GTFannotations.py --help

import sys
import pdb
import os
import subprocess
import argparse
import collections
from os.path import isfile, join
import gzip
import getopt

# typical GTF entry for one transcript of a protein coding gene
#chr1    HAVANA  gene    367640  368634  .       +       .       gene_id "ENSG00000235249.1"; transcript_id "ENSG00000235249.1"; gene_type "protein_coding"; gene_status "KNOWN"; gene_name "OR4F29"; transcript_type "protein_coding"; transcript_status "KNOWN"; transcript_name "OR4F29"; level 2; havana_gene "OTTHUMG00000002860.1";
#chr1    HAVANA  transcript      367640  368634  .       +       .       gene_id "ENSG00000235249.1"; transcript_id "ENST00000426406.1"; gene_type "protein_coding"; gene_status "KNOWN"; gene_name "OR4F29"; transcript_type "protein_coding"; transcript_status "KNOWN"; transcript_name "OR4F29-001"; level 2; tag "CCDS"; ccdsid "CCDS41220.1"; havana_gene "OTTHUMG00000002860.1"; havana_transcript "OTTHUMT00000007999.1";
#chr1    HAVANA  exon    367640  368634  .       +       .       gene_id "ENSG00000235249.1"; transcript_id "ENST00000426406.1"; gene_type "protein_coding"; gene_status "KNOWN"; gene_name "OR4F29"; transcript_type "protein_coding"; transcript_status "KNOWN"; transcript_name "OR4F29-001"; level 2; tag "CCDS"; ccdsid "CCDS41220.1"; havana_gene "OTTHUMG00000002860.1"; havana_transcript "OTTHUMT00000007999.1";
#chr1    HAVANA  CDS     367659  368594  .       +       0       gene_id "ENSG00000235249.1"; transcript_id "ENST00000426406.1"; gene_type "protein_coding"; gene_status "KNOWN"; gene_name "OR4F29"; transcript_type "protein_coding"; transcript_status "KNOWN"; transcript_name "OR4F29-001"; level 2; tag "CCDS"; ccdsid "CCDS41220.1"; havana_gene "OTTHUMG00000002860.1"; havana_transcript "OTTHUMT00000007999.1";
#chr1    HAVANA  start_codon     367659  367661  .       +       0       gene_id "ENSG00000235249.1"; transcript_id "ENST00000426406.1"; gene_type "protein_coding"; gene_status "KNOWN"; gene_name "OR4F29"; transcript_type "protein_coding"; transcript_status "KNOWN"; transcript_name "OR4F29-001"; level 2; tag "CCDS"; ccdsid "CCDS41220.1"; havana_gene "OTTHUMG00000002860.1"; havana_transcript "OTTHUMT00000007999.1";
#chr1    HAVANA  stop_codon      368595  368597  .       +       0       gene_id "ENSG00000235249.1"; transcript_id "ENST00000426406.1"; gene_type "protein_coding"; gene_status "KNOWN"; gene_name "OR4F29"; transcript_type "protein_coding"; transcript_status "KNOWN"; transcript_name "OR4F29-001"; level 2; tag "CCDS"; ccdsid "CCDS41220.1"; havana_gene "OTTHUMG00000002860.1"; havana_transcript "OTTHUMT00000007999.1";
#chr1    HAVANA  UTR     367640  367658  .       +       .       gene_id "ENSG00000235249.1"; transcript_id "ENST00000426406.1"; gene_type "protein_coding"; gene_status "KNOWN"; gene_name "OR4F29"; transcript_type "protein_coding"; transcript_status "KNOWN"; transcript_name "OR4F29-001"; level 2; tag "CCDS"; ccdsid "CCDS41220.1"; havana_gene "OTTHUMG00000002860.1"; havana_transcript "OTTHUMT00000007999.1";
#chr1    HAVANA  UTR     368595  368634  .       +       .       gene_id "ENSG00000235249.1"; transcript_id "ENST00000426406.1"; gene_type "protein_coding"; gene_status "KNOWN"; gene_name "OR4F29"; transcript_type "protein_coding"; transcript_status "KNOWN"; transcript_name "OR4F29-001"; level 2; tag "CCDS"; ccdsid "CCDS41220.1"; havana_gene "OTTHUMG00000002860.1"; havana_transcript "OTTHUMT00000007999.1";

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

def getUTRs(genedict,outname):
    print("Reading UTRs...")
    outfile=open(outname +".utrs.bed",'w')
    fiveutrlen=0
    threeutrlen=0
    for gene in genedict:
        fiveutrs=[]
        threeutrs=[]
        strand=gene.split(":")[3]
        for transcript in genedict[gene]:
            if 'start' in genedict[gene][transcript] and 'end' in genedict[gene][transcript] and 'exon_starts' in genedict[gene][transcript]:
                if strand == "+":
                    startpos = genedict[gene][transcript]['start'][0]
                    endpos = genedict[gene][transcript]['end'][1]
                if strand == "-":
                    startpos = genedict[gene][transcript]['start'][1]
                    endpos = genedict[gene][transcript]['end'][0]
                for i in range(0,len(genedict[gene][transcript]['exon_starts'])):
                    exonname=genedict[gene][transcript]['exon_starts'][i] +"-"+ genedict[gene][transcript]['exon_ends'][i]
                    exonpos=range(int(genedict[gene][transcript]['exon_starts'][i]),(int(genedict[gene][transcript]['exon_ends'][i])+1))
                    if strand == "+":
                        if int(startpos) in exonpos: fiveutrs.append(genedict[gene][transcript]['exon_starts'][i]+"-"+startpos)
                        elif int(endpos) in exonpos: threeutrs.append(endpos+"-"+genedict[gene][transcript]['exon_ends'][i])
                        elif int(startpos) > int(genedict[gene][transcript]['exon_ends'][i]): fiveutrs.append(exonname)
                        elif int(endpos) < int(genedict[gene][transcript]['exon_starts'][i]): threeutrs.append(exonname)
                    if strand == "-":
                        if int(startpos) in exonpos: fiveutrs.append(startpos+"-"+genedict[gene][transcript]['exon_ends'][i])
                        elif int(endpos) in exonpos: threeutrs.append(genedict[gene][transcript]['exon_starts'][i]+"-"+endpos)
                        elif int(startpos) < int(genedict[gene][transcript]['exon_starts'][i]): fiveutrs.append(exonname)
                        elif int(endpos) > int(genedict[gene][transcript]['exon_ends'][i]): threeutrs.append(exonname)
    fiveutrs2=list(set(fiveutrs))
    threeutrs2=list(set(threeutrs))
    fiveutrlen += len(fiveutrs2)
    threeutrlen += len(threeutrs2)
    for utr in fiveutrs2:
        outfile.write(str(gene.split(":")[1]) +"\t"+ str(utr.split("-")[0]) +"\t"+ str(utr.split("-")[1]) +"\t"+ str(gene.split(":")[0]) +"\t5UTR\t"+ str(strand) +"\n")
    for utr in threeutrs2:
        outfile.write(str(gene.split(":")[1]) +"\t"+ str(utr.split("-")[0]) +"\t"+ str(utr.split("-")[1]) +"\t"+ str(gene.split(":")[0]) +"\t3UTR\t"+ str(strand) +"\n")
    print("... " +str(fiveutrlen)+ " 5' UTRs found.")
    print("... " +str(threeutrlen)+ " 3' UTRs found.")
    outfile.close()

def getSS(genedict,outname):
    print("Reading splicesites.")
    outfile=open(outname +".splicesites.bed",'w')
    fivesslen=0
    threesslen=0
    for gene in genedict:
        fivess = []
        threess = []
        for transcript in genedict[gene]:
            for start in genedict[gene][transcript]['exon_starts']:
                if gene.split(":")[3] == "+" and start != transcript.split(":")[2].split("-")[0]: threess.append(start)
                if gene.split(":")[3] == "-" and start != transcript.split(":")[2].split("-")[0]: fivess.append(start)
            for end in genedict[gene][transcript]['exon_ends']:
                if gene.split(":")[3] == "+" and end != transcript.split(":")[2].split("-")[1]: fivess.append(end)
                if gene.split(":")[3] == "-" and end != transcript.split(":")[2].split("-")[1]: threess.append(end)
        fivess2=list(set(fivess))
        threess2=list(set(threess))
        fivesslen += len(fivess2)
        threesslen += len(threess2)
        for i in range(0,len(fivess2)):
            strand = str(gene.split(":")[3])
            outstring = str(gene.split(":")[1]) +"\t"+ str(int(fivess2[i])-3) +"\t"+ str(int(fivess2[i])+6) +"\t"+ str(gene.split(":")[0]) +"\t5ss\t"+ strand
            if strand == "-":
                outstring = str(gene.split(":")[1]) +"\t"+ str(int(fivess2[i])-7) +"\t"+ str(int(fivess2[i])+2) +"\t"+ str(gene.split(":")[0]) +"\t5ss\t"+ strand
            outfile.write(outstring +"\n")
        for i in range(0,len(threess2)):
            strand = str(gene.split(":")[3])
            outstring = str(gene.split(":")[1]) +"\t"+ str(int(threess2[i])-21) +"\t"+ str(int(threess2[i])+2) +"\t"+ str(gene.split(":")[0]) +"\t3ss\t"+ strand
            if strand == "-":
                outstring = str(gene.split(":")[1]) +"\t"+ str(int(threess2[i])-3) +"\t"+ str(int(threess2[i])+20) +"\t"+ str(gene.split(":")[0]) +"\t3ss\t"+ strand
            outfile.write(outstring +"\n")
    print("... " +str(fivesslen)+ " 5' splicesites found.")
    print("... " +str(threesslen)+ " 3' splicesites found.")
    outfile.close()

def getExons(genedict,outname):
    print("Reading exons.")
    outfile=open(outname +".exons.bed",'w')
    exonlen=0
    for gene in genedict:
        exons=[]
        strand=gene.split(":")[3]
        for transcript in genedict[gene]:
            if 'exon_starts' not in genedict[gene][transcript]: continue
            if 'start' in genedict[gene][transcript] and 'end' in genedict[gene][transcript]:
                if strand == "+":
                    startpos = genedict[gene][transcript]['start'][0]
                    endpos = genedict[gene][transcript]['end'][1]
                if strand == "-":
                    startpos = genedict[gene][transcript]['start'][1]
                    endpos = genedict[gene][transcript]['end'][0]
                for i in range(0,len(genedict[gene][transcript]['exon_starts'])):
                    exonname=genedict[gene][transcript]['exon_starts'][i] +"-"+ genedict[gene][transcript]['exon_ends'][i]
                    exonpos = range(int(genedict[gene][transcript]['exon_starts'][i]),(int(genedict[gene][transcript]['exon_ends'][i])+1))
                    if strand == "+":
                        if int(startpos) in exonpos: exons.append(startpos +"-"+ genedict[gene][transcript]['exon_ends'][i])
                        elif int(endpos) in exonpos: exons.append(genedict[gene][transcript]['exon_starts'][i] +"-"+ endpos)
                        elif int(startpos) < int(genedict[gene][transcript]['exon_starts'][i]) and int(endpos) > int(genedict[gene][transcript]['exon_ends'][i]):
                            exons.append(exonname)
                    if strand == "-":
                        if int(startpos) in exonpos: exons.append(genedict[gene][transcript]['exon_starts'][i] +"-"+ startpos)
                        elif int(endpos) in exonpos: exons.append(endpos +"-"+ genedict[gene][transcript]['exon_ends'][i])
                        elif int(startpos) > int(genedict[gene][transcript]['exon_ends'][i]) and int(endpos) < int(genedict[gene][transcript]['exon_starts'][i]):
                            exons.append(exonname)
            else:
                for i in range(0,len(genedict[gene][transcript]['exon_starts'])):
                    exons.append(genedict[gene][transcript]['exon_starts'][i] +"-"+ genedict[gene][transcript]['exon_ends'][i])
        exons2=list(set(exons))
        exonlen += len(exons2)
        for exondef in exons2:
            outfile.write(str(gene.split(":")[1]) +"\t"+ str(exondef.split("-")[0]) +"\t"+ str(exondef.split("-")[1]) +"\t"+ str(gene.split(":")[0]) +"\texon\t"+ str(strand) +"\n")
    print("... " +str(exonlen)+ " coding exons found.")
    outfile.close()

def getCEs(genedict,outname):
    print("Reading constitutive exons.")
    outfile=open(outname +".const.bed",'w')
    constlen=0
    for gene in genedict:
        consts=[]
        strand=gene.split(":")[3]
        for transcript in genedict[gene]:
            if 'exon_starts' not in genedict[gene][transcript]: continue
            exonnames=[]
            for i in range(0,len(genedict[gene][transcript]['exon_starts'])): exonnames.append(genedict[gene][transcript]['exon_starts'][i] +"-"+ genedict[gene][transcript]['exon_ends'][i])
            if not consts: consts = exonnames
            for exon in consts:
                if exon not in exonnames: consts.remove(exon)
        constlen+=len(consts)
        for exon in consts:
            outfile.write(str(strand) +"\t"+ str(exon.split("-")[0]) +"\t"+ str(exon.split("-")[1]) +"\t"+ str(gene.split(":")[0]) +"\tconst\t"+ str(strand) +"\n")
    print("... " +str(constlen)+ " constitutive exons found.")
    outfile.close()

def getIntrons(genedict,outname):
    print("Reading introns.")
    outfile=open(outname +".introns.bed",'w')
    outfiledist = open(outname +".introns_threedist.txt", 'w')
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
        for one in introns2:
            name = str(gene.split(":")[0]) +':'+ str(gene.split(":")[1]) +':'+ str(one) +':'+ str(gene.split(":")[3])
            outfile.write(str(gene.split(":")[1]) +"\t"+ str(one.split("-")[0]) +"\t"+ str(one.split("-")[1]) +"\t"+ str(name) +"\tintron\t"+ str(gene.split(":")[3]) +"\n")
            # write file with distance to 3' gene end (for kinetics scripts)
            threedist = threeend - int(one.split("-")[1])
            if strand == "-": threedist = int(one.split("-")[0]) - threeend
            outfiledist.write(str(name) +"\t"+ str(threedist) +"\n")
    print("... " +str(intronlen)+ " introns found.")
    outfile.close()
    outfiledist.close()

def getCIs(genedict,outname):
    print("Reading constitutive introns.")
    outfile=open(outname +".constintrons.bed",'w')
    outfiledist = open(outname +".constintrons_threedist.txt", 'w')
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
    parser = argparse.ArgumentParser()
    parser.add_argument('--gtf', type=str, help='gtf file from which to get annotations (full path)')
    parser.add_argument('--types',type=str, help='comma-seq list of types of annotation to output into a file (accepted: UTR,splicesites,exons,const,introns,constintrons)')
    parser.add_argument('--outname',type=str, help='basename of output file (full path)')
    args = parser.parse_args()

    genedict=getGenes(args.gtf)

#    pdb.set_trace()    

    regions=args.types
    # so far, types include: 5/3UTRs, 5/3splicesites, coding exons, constitutive exons, introns
    for type in regions.split(","):
#        print str(type)
        if type == "UTR": getUTRs(genedict,args.outname)
        if type == "splicesites": getSS(genedict,args.outname)
        if type == "exons": getExons(genedict,args.outname)
        if type == "const": getCEs(genedict,args.outname)
        if type == "introns": getIntrons(genedict,args.outname)
        if type == "constintrons": getCIs(genedict, args.outname)
