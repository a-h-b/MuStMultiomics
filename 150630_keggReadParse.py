#!/usr/bin/env python

# script to get the KEGG annotation from HMM search for each gene and sum up reads mapping to these genes by KO;
# gives for every gene the KO, together with the score of the best hit and the number of hits; 
# takes 5 inputs: - the hit files for each database (output of consolidate_hmmscan_results_justKEGG.pl);
#		  - the coverage file, a tab separated table with the number of mapped fragments per gene;
#                 optional arguments: -g, the number of genes in the dataset;
#                                     -s, the lower cut-off for the score, will be ignored it -g is provided;
#                                     -p, the cut-off for the score, if multiple hits are kept; default is only best hit;
#                                     -v, produces a list of all KOs that were detected in a file ("...validKOs.txt")
# the are 2 (or 3, if -v is set) outputs: - reads per KO
#					  - table with annotated genes, KO and score
#					  - list with assigned KOs, if -v is set
# output is written into files with the same name as the input (- the file extension), with "hits.tsv" and "hits+reads.tsv" appended

# written by Anna Heintz-Buschart (June 2015)


import os
import sys
import argparse
import math

parser = argparse.ArgumentParser(description='Select significant KOs from HMM-output.')
parser.add_argument('keggFile', help='output from consolidate_hmmscan_results_justKEGG.pl')
parser.add_argument('covFile', help='file with genes, number of reads per gene; tsv')
parser.add_argument('-g','--numberOfGenes', type=int,help='number of genes used as input to hmmer, score cut-off is calculated as log2 of this')
parser.add_argument('-s','--scoreCutoff', type=float,default=20.0,help='lower cut-off for score, defaults to 20.0, however use of -g is recommended and -s will be ignored, if -g is used')
parser.add_argument('-p','--percentageCutoff', type=float,default=80.0,help='lower cut-off for keeping multiple KOs per gene; percentage of highest scoring genes, defaults to 80.0; 100.0 would keep only best hit per gene')
parser.add_argument('-v','--validKOs', action='store_true',help='set -v, if a list of all KOs that were assigned should be written to a file')


args = parser.parse_args()
keggFile = args.keggFile
covFile = args.covFile
if args.numberOfGenes:
    sigVal = math.log(args.numberOfGenes,2)
else:
    sigVal = args.scoreCutoff
multiHitPerc = args.percentageCutoff/100

outFile = keggFile[:-3]+covFile[:-3] + "hits.tsv"
if args.validKOs:
  outFile2 = keggFile[:-3]+covFile[:-3] + "validKOs.txt"
outFile3 = keggFile[:-3]+covFile[:-3] + "hits+reads.tsv"

gene_dict = {}
kegg_file = open(keggFile, "r")
header = 1
while 1:
    linek = kegg_file.readline()
    if linek == "":
        break
    if header == 1:
        header = 0
    else:
        linek = linek.rstrip()
        tabk = linek.split("\t")
        if float(tabk[2]) >= sigVal:
            tabko, tabgene, tabscore = tabk[0].split("_")[0], tabk[1], float(tabk[2])
            if tabgene not in gene_dict:
                gene_dict[tabgene] = [[], []]
                gene_dict[tabgene][0].append(tabko)
                gene_dict[tabgene][1].append(float(tabscore))
            else:
                if tabscore >= multiHitPerc * gene_dict[tabgene][1][0]:
                    if tabscore > gene_dict[tabgene][1][0]:
                        gene_dict[tabgene][0].insert(0,tabko)
                        gene_dict[tabgene][1].insert(0,float(tabscore))
                    else:
                        gene_dict[tabgene][0].append(tabko)
                        gene_dict[tabgene][1].append(float(tabscore))
kegg_file.close()

out_file = open(outFile, "w")
out_file.write("Gene\tKO\tmaxScore\thitNumber\n")
allKOs = []
gene_dict_tidy = {}
for item in gene_dict:
    gene = item
    priKOs = []
    hN = 0
    score = gene_dict[item][1][0]
    for KOind in range(len(gene_dict[item][0])):
        if gene_dict[item][1][KOind] >= multiHitPerc * score and gene_dict[item][0][KOind] not in priKOs:
            priKOs.append(gene_dict[item][0][KOind])
            if gene_dict[item][0][KOind] not in allKOs:
                allKOs.append(gene_dict[item][0][KOind])
    koIDs = ";".join(priKOs)
    hN = len(priKOs)
    gene_dict_tidy[gene] = priKOs
    out_file.write(gene + "\t" + koIDs + "\t" + str(score) + "\t" + str(hN) + "\n")
out_file.close()
if args.validKOs:
    out_file2 = open(outFile2, "w")
    for kos in allKOs:
        out_file2.write(str(kos)+"\n")
    out_file2.close()
print(len(allKOs))

cov_file = open(covFile, "r")
allKOVals = [allKOs,[0.0]*len(allKOs)] #0:name, 1:reads 
otherKOs = [] #reads
header = 1
while 1:
    linec = cov_file.readline()
    if linec == "" and header==0:
        break
    if header == 1:
        header = 0
    else:
        linec = linec.rstrip()
	tabc = linec.split("\t") #0: name, 1:reads
	gene = tabc[0]
	reads = float(tabc[1])
	if gene in gene_dict_tidy:
	    funLen = len(gene_dict_tidy[gene])
	    for KO in gene_dict_tidy[gene]:
		indx = allKOVals[0].index(KO)
		allKOVals[1][indx] += reads/funLen
	else:
	    otherKOs.append(reads)
out_file3 = open(outFile3,"w")
out_file3.write("KO\treads\n")
for i in range(len(allKOVals[0])):
    if allKOVals[1][i] > 0:
	out_file3.write(allKOVals[0][i]+"\t"+str(allKOVals[1][i])+"\n")
out_file3.write("other\t"+str(sum(otherKOs))+"\n")        
    
