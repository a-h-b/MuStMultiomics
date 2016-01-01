#!/usr/bin/env python


# script to get the best annotation from HMM search with 5 datbases (KEGG, metaCyc, Swiss-Prot, Pfam and TIGR Pfam) for each gene
#         and sum up reads mapping to these genes by function;
# takes 7 inputs: - the hit files for each database (output of consolidate_hmmscan_results.pl) in this order: KEGG, metaCyc, Swiss-Prot, Pfam, TIGR Pfam;
#		  - the coverage file, a tab separated table with the number of mapped fragments per gene;
#                 -g, the number of genes in the dataset (non-optional in this script);
# the is only 1 output, a table with reads per function
# output is written into a file with the same name as the input coverage file (- the file extension) with "besthitsAllDB" prepended and "reads.tsv" appended

# written by Anna Heintz-Buschart (March 2015)

 
import os
import sys
import argparse
import math

parser = argparse.ArgumentParser(description='Select significant KOs from HMM-output.')
parser.add_argument('koFile', help='KEGG output files from consolidate_hmmscan_results.pl')
parser.add_argument('mcFile', help='metaCyc output files from consolidate_hmmscan_results.pl')
parser.add_argument('spFile', help='Swiss-Prot output files from consolidate_hmmscan_results.pl')
parser.add_argument('pfFile', help='Pfam output files from consolidate_hmmscan_results.pl')
parser.add_argument('tiFile', help='TIGR Pfam output files from consolidate_hmmscan_results.pl')
parser.add_argument('covFile', help='file with genes, number of reads; tsv')
parser.add_argument('-g','--numberOfGenes', type=int,help='number of genes used as input to hmmer, score cut-off is calculated as log2 of this')


args = parser.parse_args()
koFile = args.koFile
mcFile = args.mcFile
spFile = args.spFile
pfFile = args.pfFile
tiFile = args.tiFile
covFile = args.covFile
if args.numberOfGenes:
    sigVal = math.log(args.numberOfGenes,2)
annN = "ID"

outFile = "besthitsAllDB."+covFile[:-7] + "reads.tsv"

gene_dict = {}
hmm_file = open(koFile, "r")
header = 1
while 1:
    linek = hmm_file.readline()
    if linek == "":
        break
    if header == 1:
        header = 0
    else:
        linek = linek.rstrip()
        tabi = linek.split("\t")
        if float(tabi[2]) >= sigVal:
            tabid, tabgene, tabscore = "KEGG:"+tabi[0].split("_")[0], tabi[1], float(tabi[2])
            if tabgene not in gene_dict:
                gene_dict[tabgene] = [[], []]
                gene_dict[tabgene][0].append(tabid)
                gene_dict[tabgene][1].append(float(tabscore))
            else:
                if tabscore >= gene_dict[tabgene][1][0]:
                    if tabscore > gene_dict[tabgene][1][0]:
                        gene_dict[tabgene][0].insert(0,tabid)
                        gene_dict[tabgene][1].insert(0,float(tabscore))
                    else:
                        gene_dict[tabgene][0].append(tabid)
                        gene_dict[tabgene][1].append(float(tabscore))
hmm_file.close()

hmm_file = open(mcFile, "r")
header = 1
while 1:
    linek = hmm_file.readline()
    if linek == "":
        break
    if header == 1:
        header = 0
    else:
        linek = linek.rstrip()
        tabi = linek.split("\t")
        if float(tabi[2]) >= sigVal:
            tabid, tabgene, tabscore = "metaCyc:"+tabi[0].split("_")[0], tabi[1], float(tabi[2])
            if tabgene not in gene_dict:
                gene_dict[tabgene] = [[], []]
                gene_dict[tabgene][0].append(tabid)
                gene_dict[tabgene][1].append(float(tabscore))
            else:
                if tabscore >= gene_dict[tabgene][1][0]:
                    if tabscore > gene_dict[tabgene][1][0]:
                        gene_dict[tabgene][0].insert(0,tabid)
                        gene_dict[tabgene][1].insert(0,float(tabscore))
                    else:
                        gene_dict[tabgene][0].append(tabid)
                        gene_dict[tabgene][1].append(float(tabscore))
hmm_file.close()

hmm_file = open(spFile, "r")
header = 1
while 1:
    linek = hmm_file.readline()
    if linek == "":
        break
    if header == 1:
        header = 0
    else:
        linek = linek.rstrip()
        tabi = linek.split("\t")
        if float(tabi[2]) >= sigVal:
            tabid, tabgene, tabscore = "swissProt:"+tabi[0].split("_")[0], tabi[1], float(tabi[2])
            if tabgene not in gene_dict:
                gene_dict[tabgene] = [[], []]
                gene_dict[tabgene][0].append(tabid)
                gene_dict[tabgene][1].append(float(tabscore))
            else:
                if tabscore >= gene_dict[tabgene][1][0]:
                    if tabscore > gene_dict[tabgene][1][0]:
                        gene_dict[tabgene][0].insert(0,tabid)
                        gene_dict[tabgene][1].insert(0,float(tabscore))
                    else:
                        gene_dict[tabgene][0].append(tabid)
                        gene_dict[tabgene][1].append(float(tabscore))
hmm_file.close()

hmm_file = open(pfFile, "r")
header = 1
while 1:
    linek = hmm_file.readline()
    if linek == "":
        break
    if header == 1:
        header = 0
    else:
        linek = linek.rstrip()
        tabi = linek.split("\t")
        if float(tabi[2]) >= sigVal:
            tabid, tabgene, tabscore = "Pfam:"+tabi[0], tabi[1], float(tabi[2])
            if tabgene not in gene_dict:
                gene_dict[tabgene] = [[], []]
                gene_dict[tabgene][0].append(tabid)
                gene_dict[tabgene][1].append(float(tabscore))
            else:
                if tabscore >= gene_dict[tabgene][1][0]:
                    if tabscore > gene_dict[tabgene][1][0]:
                        gene_dict[tabgene][0].insert(0,tabid)
                        gene_dict[tabgene][1].insert(0,float(tabscore))
                    else:
                        gene_dict[tabgene][0].append(tabid)
                        gene_dict[tabgene][1].append(float(tabscore))
hmm_file.close()

hmm_file = open(tiFile, "r")
header = 1
while 1:
    linek = hmm_file.readline()
    if linek == "":
        break
    if header == 1:
        header = 0
    else:
        linek = linek.rstrip()
        tabi = linek.split("\t")
        if float(tabi[2]) >= sigVal:
            tabid, tabgene, tabscore = "TIGR:"+tabi[0].split("_")[0], tabi[1], float(tabi[2])
            if tabgene not in gene_dict:
                gene_dict[tabgene] = [[], []]
                gene_dict[tabgene][0].append(tabid)
                gene_dict[tabgene][1].append(float(tabscore))
            else:
                if tabscore >= gene_dict[tabgene][1][0]:
                    if tabscore > gene_dict[tabgene][1][0]:
                        gene_dict[tabgene][0].insert(0,tabid)
                        gene_dict[tabgene][1].insert(0,float(tabscore))
                    else:
                        gene_dict[tabgene][0].append(tabid)
                        gene_dict[tabgene][1].append(float(tabscore))
hmm_file.close()


allIDs = []
gene_dict_tidy = {}
for item in gene_dict:
    gene = item
    priIDs = []
    hN = 0
    score = gene_dict[item][1][0]
    for IDind in range(len(gene_dict[item][0])):
        if gene_dict[item][1][IDind] >= score and gene_dict[item][0][IDind] not in priIDs:
            priIDs.append(gene_dict[item][0][IDind])
            if gene_dict[item][0][IDind] not in allIDs:
                allIDs.append(gene_dict[item][0][IDind])
    koIDs = ";".join(priIDs)
    hN = len(priIDs)
    gene_dict_tidy[gene] = priIDs
print(len(allIDs))

cov_file = open(covFile, "r")
allIDVals = [allIDs,[0.0]*len(allIDs)] #0:name, 1:reads
otherIDs = [] #reads
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
	    for ID in gene_dict_tidy[gene]:
		indx = allIDVals[0].index(ID)
		allIDVals[1][indx] += reads/funLen
	else:
	    otherIDs.append(reads)
out_file = open(outFile,"w")
out_file.write("DB:ID\treads\n")
for i in range(len(allIDVals[0])):
    if allIDVals[1][i] > 0:
	out_file.write(allIDVals[0][i]+"\t"+str(allIDVals[1][i])+"\n")
out_file.write("other\t"+str(sum(otherIDs))+"\n")
