#!/usr/bin/env python

# script to get the KEGG annotation from HMM search for each gene, limited to KOs that are part of a reconstructed network;
# gives for every gene the node, together with the score of the best hit and the number of hits
# takes 4 inputs: the hit files for each database (output of consolidate_hmmscan_results_justKEGG.pl); 
#                 optional arguments: -g, the number of genes in the dataset;
#                                     -s, the lower cut-off for the score, will be ignored it -g is provided;
#                                     -p, the cut-off for the score, if multiple hits are kept; default is only best hit;
#                                     -v, produces a list of all KOs that were detected in a file ("...validKOs.txt")
# a file with all nodes and KOs in the network must be provided (1st column: node, 2nd column: KO); path is hard coded below
# output is written into a file with the same name as the input (- the file extension), with "hitsNodes.tsv" appended
# written by Anna Heintz-Buschart (July 2015)

import os
import sys
import argparse
import math

nodesFile = "150705_KOs_in_NW.tsv"

parser = argparse.ArgumentParser(description='Select significant KOs from HMM-output.')
parser.add_argument('keggFile', help='output from consolidate_hmmscan_results_justKEGG.pl')
parser.add_argument('-g','--numberOfGenes', type=int,help='number of genes used as input to hmmer, score cut-off is calculated as log2 of this')
parser.add_argument('-s','--scoreCutoff', type=float,default=20.0,help='lower cut-off for score, defaults to 20.0, however use of -g is recommended and -s will be ignored, if -g is used')
parser.add_argument('-p','--percentageCutoff', type=float,default=100.0,help='lower cut-off for keeping multiple KOs per gene; percentage of highest scoring genes, defaults to 100.0; 100.0 would keep only best hit per gene')
parser.add_argument('-v','--validKOs', action='store_true',help='set -v, if a list of all KOs that were assigned should be written to a file')

args = parser.parse_args()
keggFile = args.keggFile
if args.numberOfGenes:
    sigVal = math.log(args.numberOfGenes,2)
else:
    sigVal = args.scoreCutoff
multiHitPerc = args.percentageCutoff/100

outFile = keggFile[:-3] + "hitsNodes.tsv"
if args.validKOs:
    outFile2 = keggFile[:-3] + "validKOs.txt"

nodes_dict = {}
nodes_file = open(nodesFile,"r")
header = 1
while 1:
    linek = nodes_file.readline()
    if linek == "":
        break
    if header == 1:
        header = 0
    else:
        linek = linek.rstrip()
        tabk = linek.split("\t")
        nodes_dict[tabk[1]] = tabk[0]
nodes_file.close()

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
            if tabko in nodes_dict:
                if tabgene not in gene_dict:
                    gene_dict[tabgene] = [[], []]
                    gene_dict[tabgene][0].append(nodes_dict[tabko])
                    gene_dict[tabgene][1].append(float(tabscore))
                else:
                    if tabscore >= multiHitPerc * gene_dict[tabgene][1][0]:
                        if tabscore > gene_dict[tabgene][1][0]:
                            gene_dict[tabgene][0].insert(0,nodes_dict[tabko])
                            gene_dict[tabgene][1].insert(0,float(tabscore))
                        else:
                            gene_dict[tabgene][0].append(nodes_dict[tabko])
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
