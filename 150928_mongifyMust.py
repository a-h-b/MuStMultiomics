#!/usr/bin/env python


#this script is not meant to run anywhere but the environment it was originally written for, as it is full of hard-coded paths.
# consider it documentation :-)
# written by Anna Heintz-Buschart, this version is the final version used in the MuSt project, September 2015

import os
import sys
sys.path.append('/home/users/aheintzbuschart/lib/')
import argparse
import re
import numpy
import collections
from pymongo import MongoClient

# get arguments from call
parser = argparse.ArgumentParser(description='Feed mongo DB with contig and gene annotations')
parser.add_argument('-f','--family', type=str,help='family ID ie must_m_01')
parser.add_argument('-l','--lib', type=str,help='sample ID ie M1.1-V1')

args = parser.parse_args()
FAM = args.family #family membership of this sample
LIB = args.lib # combined assembly / metaT ID of this sample

#function to be used to simplify lists of taxa
def enumerateTaxa(taxList):
    kp = collections.Counter(taxList).most_common()
    init = 1
    taxString = ""
    for tax,occ in kp:
        if init == 1:
            init = 0
            taxString = tax + "(" + str(occ) + ");"
        else:
            taxString += tax + "(" + str(occ) + ");"
    taxString = taxString.rstrip(";")
    return taxString


# read file linking combined assembly / metaT ID to metaG ID to get metaG ID
idFile = "/work/projects/ecosystem_biology/MUST/CombinedAssembly/Annotation/sequenceInfo/" + FAM + "/ids"
print "Reading ids from ", idFile 
id_file = open(idFile, "r")
while 1:
    line = id_file.readline()
    if line == "":
        break
    else:
        line = line.rstrip()
        tab = line.split("\t") # 0: DNAid 1: combi ID
        if tab[1] == LIB:
            dLIB = tab[0]
id_file.close()
print "DNA library for this sample is ", dLIB

# read file containing all samples belonging the individuals in this family and get names of metaG ids for the other samples of the same individual
visit_list = []
visitFile = "/work/projects/ecosystem_biology/MUST/CombinedAssembly/Annotation/mapping/" + FAM + "/visits"
print "Reading other sample IDs from ", visitFile
visit_file = open(visitFile, "r")
while 1:
    line = visit_file.readline()
    if line == "":
        break
    else:
        line = line.rstrip()
        if dLIB in line:
            visit_list = line.split(" ")
            for vis in visit_list:
                if dLIB == vis:
                    visit_list.remove(vis)
visit_file.close()
print "other DNA libraries for this individual are ", ",".join(visit_list)

# read file linking proteomics results to combined ID
proteomeIdFile = "/work/projects/ecosystem_biology/MUST/MUST-metaP/GIGAsearches/" + FAM + "/GIGA_ind_combined_ids"
print "Reading ids and file names for proteomics from ", proteomeIdFile
proteomeId_file = open(proteomeIdFile, "r")
while 1:
    line = proteomeId_file.readline()
    if line == "":
        break
    else:
        line = line.rstrip()
        tab = line.split(",") # 0: combiID, 1: individual, 2:proteomics prefix
        if tab[0] == LIB:
            pLIB = tab[2]
pFile = "/work/projects/ecosystem_biology/MUST/MUST-metaP/GIGAsearches/exportedResultsQuan/" + pLIB + "_TargetProtein.txt"
pgFile = "/work/projects/ecosystem_biology/MUST/MUST-metaP/GIGAsearches/exportedResultsQuan/" + pLIB + "_TargetProteinGroup.txt"
visit_file.close()
print "proteomics results for this individual are in files ", pFile , " and ", pgFile


# read file of contig information and make dictionary with all contigs
contig_dict = {}
header = 1
contigFile = "/work/projects/ecosystem_biology/MUST/CombinedAssembly/Annotation/sequenceInfo/" + FAM + "/" + LIB + "/contigs.length.tsv"
print "Reading contig information from ", contigFile
contig_file = open(contigFile, "r")
while 1:
    line = contig_file.readline()
    if line == "":
        break
    if header == 1:
        header = 0
    else:
        line = line.rstrip()
        tab = line.split("\t") # 0:sequenceID, 1:length, 2:GCperc
        contig_dict[tab[0]] = {'sample': LIB,'length' : int(tab[1]), 'GCperc' : float(tab[2]), 'aveCov' : 0, 'varPerMB' : 0, 'cluster' : "S"}
contig_file.close()

# read file of DNA coverage of the same sample for contigs and update contig dictionary
header = 1
covDNAcontigFile = "/work/projects/ecosystem_biology/MUST/CombinedAssembly/Annotation/mapping/" + FAM + "/" + LIB + "/" + dLIB + ".DNAonContigs.cov.tsv"
print "Reading coverage data for DNA (this sample) on contigs from ", covDNAcontigFile
covDNAcontig_file = open(covDNAcontigFile, "r")
while 1:
    line = covDNAcontig_file.readline()
    if line == "":
        break
    if header == 1:
        header = 0
    else:
        line = line.rstrip()
        tab = line.split("\t") # 0:SequenceID,	1:Reference.length,	2: Average.coverage,	3:Covered.length
        contig_dict[tab[0]]['aveCov'] = float(tab[2])
covDNAcontig_file.close()

# read  DNA coverages of the other samples of the same individual for contigs and update contig dictionary
for other in visit_list:
    visit = "aveCovV" + other[other.index('V')+1]
    header = 1
    covDNAcontigFile = "/work/projects/ecosystem_biology/MUST/CombinedAssembly/Annotation/mapping/" + FAM + "/" + LIB + "/" + other + ".DNAonContigs.cov.tsv"
    print "Reading coverage data for DNA of ", other," on contigs from ", covDNAcontigFile
    covDNAcontig_file = open(covDNAcontigFile, "r")
    while 1:
        line = covDNAcontig_file.readline()
        if line == "":
            break
        if header == 1:
            header = 0
        else:
            line = line.rstrip()
            tab = line.split("\t") # 0:SequenceID,	1:Reference.length,	2: Average.coverage,	3:Covered.length
            contig_dict[tab[0]][visit] = float(tab[2])
    covDNAcontig_file.close()

# read variants file
header = 1
varContig_dict = {} 
varGene_dict = {}
varFile = "/work/projects/ecosystem_biology/MUST/CombinedAssembly/Annotation/variants/" + FAM + "/" + LIB + "/variantsAnnotation.tsv"
print "Reading variant information from ", varFile
var_file = open(varFile, "r")
while 1:
    line = var_file.readline()
    if line == "":
        break
    if header == 1:
        header = 0
    else:
        line = line.rstrip()
        tab = line.split("\t") # 0:contig, 1:pos, 2:type, 3:diff, 4:gene, 5:frameshift, 6:totalCov, 7:totalVarCov, 8:DNACov, 9:DNAvarCov
        ### 10:RNACov, 11:RNAvarCov
        if float(tab[8]) > 0:
            varRelColv = float(tab[9])/float(tab[8])
        else:
            varRelColv = 0
        if tab[0] not in varContig_dict:
            varContig_dict[tab[0]] = {'varPerMB' : 1.0, 'varRelCov' : [varRelColv], 'varPos' : [int(tab[1])]}
        else:
            varContig_dict[tab[0]]['varPerMB'] += 1
            varContig_dict[tab[0]]['varRelCov'].append(varRelColv)
            varContig_dict[tab[0]]['varPos'].append(int(tab[1]))
        # the positions come as a list now
        # the relative variant coverage is different from R but more correct now
        if tab[4] not in varGene_dict:
            varGene_dict[tab[4]] = 1.0
        else:
            varGene_dict[tab[4]] += 1
var_file.close()

# update contig dictionary with variant data
for var in varContig_dict:
    contig_dict[var]['varPerMB'] = varContig_dict[var]['varPerMB'] *1000000 / contig_dict[var]['length']
    contig_dict[var]['varRelCovMean'] = numpy.mean(varContig_dict[var]['varRelCov'])
    contig_dict[var]['varPos'] = varContig_dict[var]['varPos']
varContig_dict = {}            

# read Kraken data and update contig dictionary
header = 1
krakenFile = "/work/projects/ecosystem_biology/MUST/CombinedAssembly/Annotation/kraken/" + FAM + "/" + LIB + "/" + LIB + ".contigs.kraken_BacVirGenBank40GBannoUnambig.tsv"
print "Reading contig Kraken information from ", krakenFile
kraken_file = open(krakenFile, "r")
while 1:
    line = kraken_file.readline()
    if line == "":
        break
    if header == 1:
        header = 0
    else:
        line = line.rstrip()
        tab = line.split("\t") # 0:contig, 1:length, 2:divSpec, 3:divGen, 4:divFam, 5:divOrd, 6:divClass, 7:divPhylum, 8:divKingdom,
        ### 9:annotationLevel, 10:species, 11:genus, 12:family, 13:order, 14:class, 15:phylum, 16:kingdom
        contig_dict[tab[0]]['krakenAnnotationLevel'] = tab[9]
        if tab[10] != "unknown":
            contig_dict[tab[0]]['krakenSpecies'] = tab[10]
        if tab[11] != "unknown":
            contig_dict[tab[0]]['krakenGenus'] = tab[11]
        if tab[12] != "unknown":
            contig_dict[tab[0]]['krakenFamily'] = tab[12]
        if tab[13] != "unknown":
            contig_dict[tab[0]]['krakenOrder'] = tab[13]
        if tab[14] != "unknown":
            contig_dict[tab[0]]['krakenClass'] = tab[14]
        if tab[15] != "unknown":
            contig_dict[tab[0]]['krakenPhylum'] = tab[15]
        if tab[16] != "unknown":
            contig_dict[tab[0]]['krakenKingdom'] = tab[16]
        # I only use annotated contigs; the names of the fields are different from R, because I cannot have "." in the keys
kraken_file.close()

# read contig names for coordinates
coord_names = []
coordNameFile = "/work/projects/ecosystem_biology/MUST/CombinedAssembly/Annotation/Bmaps/" + FAM + "/" + LIB + "/contigs.1000.rRNAcut.names.txt"
print "Reading contig names for BHSNE coordinates from ", coordNameFile
coordName_file = open(coordNameFile, "r")
while 1:
    line = coordName_file.readline()
    if line == "":
        break
    else:
        line = line.rstrip() # names
        line = line.replace(">","")
        coord_names.append(line)
coordName_file.close()

# read coordinates
i = 0
coordFile = "/work/projects/ecosystem_biology/MUST/CombinedAssembly/Annotation/Bmaps/" + FAM + "/" + LIB + "/contigs.1000.rRNAcut.fa_5mer_clr.coords"
print "Reading contig BHSNE coordinates from ", coordFile
coord_file = open(coordFile, "r")
while 1:
    line = coord_file.readline()
    if line == "":
        break
    else:
        line = line.rstrip()
        tab = line.split(",") # 0:x, 1:y
        contig_dict[coord_names[i]]['coords'] = [float(tab[0]),float(tab[1])] #different from R as both values come as one list now
        i += 1
coord_file.close()

#read cluster membership
header = 1
clusterFile = "/work/projects/ecosystem_biology/MUST/CombinedAssembly/Annotation/Bmaps/" + FAM + "/" + LIB + "/contigs2clusters.tsv"
print "Reading contig cluster membership from ", clusterFile
cluster_file = open(clusterFile, "r")
while 1:
    line = cluster_file.readline()
    if line == "":
        break
    if header == 1:
        header = 0
    else:
        line = line.rstrip()
        tab = line.split("\t") # 0:contig, 1:class (contig membership)
        contig_dict[tab[0]]['cluster'] = tab[1]
cluster_file.close()
contigNum = len(contig_dict)
print "gathered information on ", contigNum, "contigs"

# read file of protein predictions and make dictionary with all proteins
protein_dict = {}
protFile = "/work/projects/ecosystem_biology/MUST/CombinedAssembly/Annotation/sequenceInfo/" + FAM + "/" + LIB + "/gene.prediction.assembly.Prodigal.500/contig.Prodigal.tab"
print "Reading protein information from ", protFile
prot_file = open(protFile, "r")
while 1:
    line = prot_file.readline()
    if line == "":
        break
    else:
        line = line.rstrip()
        tab = line.split("\t") # 0:gene, 1:sense, 2:start, 3:end, 4:length, 5:startCodon, 6:stopCodon, 7:completeness
        protein_dict[tab[0]] = {'sense' : tab[1], 'start' : int(tab[2]), 'end' : int(tab[3]), 'length' : int(tab[4]), 'startCodon' : tab[5], 'stopCodon' : tab[6], 'completeness' : tab[7], 'kind' : "protein", 'aveCovDNA' : 0, 'aveCovRNAfw' : 0, 'readsRNAfw' : 0, 'aveCovRNArc' : 0, 'varPerMB' : 0}
prot_file.close()

# read file of RNA predictions and make dictionary with all RNAs
rna_dict = {}
header = 1
rnaFile = "/work/projects/ecosystem_biology/MUST/CombinedAssembly/Annotation/hmms/" + FAM + "/" + LIB + "/rRNAgenes.bac.tab"
print "Reading RNA information from ", rnaFile
rna_file = open(rnaFile, "r")
while 1:
    line = rna_file.readline()
    if line == "":
        break
    if header == 1:
        header = 0
    else:
        line = line.rstrip()
        tab = line.split("\t") # 0:contig, 1:gene, 2:sense, 3:length, 4:start, 5:end, 6:completeness, 7:kind
        rna_dict[tab[1]] = {'sense' : tab[2], 'start' : int(tab[4]), 'end' : int(tab[5]), 'length' : int(tab[3]), 'completeness' : tab[6], 'kind' : tab[7], 'aveCovDNA' : 0, 'aveCovRNAfw' : 0, 'readsRNAfw' : 0, 'aveCovRNArc' : 0, 'varPerMB' : 0}
rna_file.close()

# read file of DNA coverage for genes and update protein or RNA dictionary
header = 1
covDNAgeneFile = "/work/projects/ecosystem_biology/MUST/CombinedAssembly/Annotation/mapping/" + FAM + "/" + LIB + "/" + dLIB + ".DNAonGenesrRNA.cov.tsv"
print "Reading coverage data for DNA on genes from ", covDNAgeneFile
covDNAgene_file = open(covDNAgeneFile, "r")
while 1:
    line = covDNAgene_file.readline()
    if line == "":
        break
    if header == 1:
        header = 0
    else:
        line = line.rstrip()
        tab = line.split("\t") # 0:SequenceID,	1:Reference.length,	2: Average.coverage,	3:Covered.length
        if tab[0] in protein_dict:
            protein_dict[tab[0]]['aveCovDNA'] = float(tab[2])
        elif tab[0] in rna_dict:
            rna_dict[tab[0]]['aveCovDNA'] = float(tab[2])
covDNAgene_file.close()

# read file of RNA forward coverage for genes and update protein or RNA dictionary
header = 1
covRNAgeneFile = "/work/projects/ecosystem_biology/MUST/CombinedAssembly/Annotation/mapping/" + FAM + "/" + LIB + "/" + LIB + ".RNAonGenesrRNA.cov.tsv"
print "Reading forward coverage data for RNA on genes from ", covRNAgeneFile
covRNAgene_file = open(covRNAgeneFile, "r")
while 1:
    line = covRNAgene_file.readline()
    if line == "":
        break
    if header == 1:
        header = 0
    else:
        line = line.rstrip()
        tab = line.split("\t") # 0:SequenceID,	1:Reference.length,	2: Average.coverage,	3:Covered.length
        if tab[0] in protein_dict:
            protein_dict[tab[0]]['aveCovRNAfw'] = float(tab[2])
        elif tab[0] in rna_dict:
            rna_dict[tab[0]]['aveCovRNAfw'] = float(tab[2])
covRNAgene_file.close()

# read file of RNA forward mapping reads for genes and update protein or RNA dictionary
header = 1
mapRNAgeneFile = "/work/projects/ecosystem_biology/MUST/CombinedAssembly/Annotation/mapping/" + FAM + "/" + LIB + "/" + LIB + ".RNAonGenesrRNA.mapped.tsv"
print "Reading forward mapping reads numbers for RNA on genes from ", mapRNAgeneFile
mapRNAgene_file = open(mapRNAgeneFile, "r")
while 1:
    line = mapRNAgene_file.readline()
    if line == "":
        break
    if header == 1:
        header = 0
    else:
        line = line.rstrip()
        tab = line.split("\t") # 0:SequenceID,	1:Reads
        if tab[0] in protein_dict:
            protein_dict[tab[0]]['readsRNAfw'] = int(tab[1])
        elif tab[0] in rna_dict:
            rna_dict[tab[0]]['readsRNAfw'] = int(tab[1])
mapRNAgene_file.close()


# read file of RNA reverse complement coverage for genes and update protein or RNA dictionary
header = 1
covRNArcgeneFile = "/work/projects/ecosystem_biology/MUST/CombinedAssembly/Annotation/mapping/" + FAM + "/" + LIB + "/" + LIB + ".revcompRNAonGenesrRNA.cov.tsv"
print "Reading reverse complement coverage data for RNA on genes from ", covRNArcgeneFile
covRNArcgene_file = open(covRNArcgeneFile, "r")
while 1:
    line = covRNArcgene_file.readline()
    if line == "":
        break
    if header == 1:
        header = 0
    else:
        line = line.rstrip()
        tab = line.split("\t") # 0:SequenceID,	1:Reference.length,	2: Average.coverage,	3:Covered.length
        if tab[0] in protein_dict:
            protein_dict[tab[0]]['aveCovRNArc'] = float(tab[2])
        elif tab[0] in rna_dict:
            rna_dict[tab[0]]['aveCovRNArc'] = float(tab[2])
covRNArcgene_file.close()

# update protein or RNA dictionary with variant data
for var in varGene_dict:
    if var in protein_dict:
            protein_dict[var]['varPerMB'] = varGene_dict[var] *1000000 / protein_dict[var]['length']
    elif var in rna_dict:
            rna_dict[var]['varPerMB'] = varGene_dict[var] *1000000 / rna_dict[var]['length']
varGene_dict = {}            

# read protein group information
header = 1
cntPG = 0
pg_dict = {}
print "Reading protein group information from ", pgFile
pg_file = open(pgFile, "r")
while 1:
    line = pg_file.readline()
    if line == "":
        break
    if header == 1:
        header = 0
    else:
        line = line.replace('"', '').rstrip()
        tab = line.split("\t") # 0:Protein Group ID, 1:# Proteins, 2: # Unique Peptides, 3:# Peptides, 4:# PSMs, 5:Group Description
                               ### 6:Areas:
        if int(tab[2]) > 0:
            pg_dict[tab[0]] = {'membersNum':int(tab[1])}
            cntPG += 1
pg_file.close()
print "Read protein group information for ", cntPG, " protein groups with >= 1 unique peptide"

# read protein information
header = 1
cntP = 0
p_dict = {}
print "Reading protein information from ", pFile
p_file = open(pFile, "r")
while 1:
    line = p_file.readline()
    if line == "":
        break
    if header == 1:
        header = 0
    else:
        line = line.replace('"', '').rstrip()
        tab = line.split("\t") # 0:Master, 1:Unique Sequence ID, 2:Protein Group IDs, 3:Accession, 4:Description, 5:Coverage,
                               ### 6:# Peptides, 7:# PSMs, 8:# Protein Unique Peptides, 9:# Unique Peptides, 10:# Protein Groups, 11:# AAs
                               ### 12: MW [kDa], 13: calc. pI, 14: Modifications, 15: Areas: F2: Sample, 16: emPAI, 17: Score Sequest HT
                               ### 18: Coverage Sequest HT, 19: # Peptides Sequest HT, 20: # PSMs Sequest HT
        if int(tab[3]) < 30000000: #this means the protein is not identified from the human database
            cntP += 1
            pgs = tab[2].split(";")
            pgList = []
            for pg in pgs:
                if pg in pg_dict:
                    pgList.append(pg)
            if len(pgList) != 0:
                if tab[15] == "":
                    tab[15] = 0
                p_dict[tab[3]] = {'desc' : tab[4], 'pgID' : pgList, 'area': float(tab[15]), 'peptides' : int(tab[6]), 'coverage' : float(tab[5])}
p_file.close()
print "Read protein information for ", cntP, " proteins"
print "included ", len(p_dict), " in dictionary."

# update protein dictionary with proteomics results
varProt_dict = {}
delList = []
protCnt = 0
protCntb = 0
pCnt = 0
prots = []
for p in p_dict:
    pCnt += 1
    idAs = []
    for pgm in p_dict[p]['pgID']:
        if pg_dict[pgm]['membersNum'] == 1: # uniquely identified microbial proteins
            protName = p_dict[p]['desc']
            if "variant" in protName:
                idAs += [protName.split(".")[1]]
                protName = protName.split(".")[0]
            else:
                idAs += ["ref"]
            if protName in protein_dict:
                prots.append(protName)
                if 'proteinIdentification' not in protein_dict[protName]:
                    protein_dict[protName]['proteinIdentification'] = "uniquely"
                    protein_dict[protName]['proteinIdentificationAs'] = [idAs]
                    protein_dict[protName]['proteinArea'] = [p_dict[p]['area']]
                    protein_dict[protName]['peptides'] = [p_dict[p]['peptides']]
                    protein_dict[protName]['proteinCoverage'] = [p_dict[p]['coverage']]
                    protCnt += 1
                    protCntb += 1
                else:
                    protCntb += 1
                    protein_dict[protName]['proteinIdentificationAs'] += [idAs]
                    protein_dict[protName]['proteinArea'] += [p_dict[p]['area']]
                    protein_dict[protName]['peptides'] += [p_dict[p]['peptides']]
                    protein_dict[protName]['proteinCoverage'] += [p_dict[p]['coverage']]
            delList.append(pgm)
        else: # put putatively identified proteins into protein group dictionary
            if 'allMembers' not in pg_dict[pgm]:
                pg_dict[pgm]['allMembers'] = [p]
            else:
                pg_dict[pgm]['allMembers'] += [p]
print "integrated ", protCnt, " proteins making up protein groups"
print "processed proteins ", pCnt 
print "protein groups with one microbial protein: ", len(delList)
delList = []
for pg in pg_dict: #remove all proteins with a single member (left over human)
    if pg_dict[pg]['membersNum'] == 1:
        delList.append(pg)
print "removed protein groups with one protein: ", len(delList)
for pg in delList:
    del pg_dict[pg]
print "protein groups with more than one member: ", len(pg_dict)

huoCnt = 0
pgCnt = 0
pgCntb = 0
pgCnt2 = 0
for pg in pg_dict: # find protein groups that represent only one protein with variants
    potVar = []
    idAs = []
    areas = []
    peptides = []
    protCov = []
    if 'allMembers' in pg_dict[pg]: #do only for those protein groups that contain at least one microbial protein
        for p in pg_dict[pg]['allMembers']:
            protName = ""
            pn = p_dict[p]['desc']
            potVar.append(re.sub(r'.variant[0-9]+$',"",pn))
            if "variant" not in pn:
                if "ref" not in idAs:
                    idAs += ["ref"]
            else:
                idAs += [pn.split(".")[1]]
        if len(set(potVar)) == 1:
            pgCntb += 1
            protName = list(set(potVar))[0]
            for p in pg_dict[pg]['allMembers']:
                areas.append(p_dict[p]['area'])
                peptides.append(p_dict[p]['peptides'])
                protCov.append(p_dict[p]['coverage'])
            if protName in protein_dict:
                prots.append(protName)
                if 'proteinIdentification' not in protein_dict[protName]:
                    if len(potVar) == pg_dict[pg]['membersNum']:
                        protein_dict[protName]['proteinIdentification'] = "uniquely"
                    else:
                        protein_dict[protName]['proteinIdentification'] = "putatively"
                    protein_dict[protName]['proteinIdentificationAs'] = [idAs]
                    protein_dict[protName]['proteinArea'] = areas
                    protein_dict[protName]['peptides'] = peptides
                    protein_dict[protName]['proteinCoverage'] = protCov
                    protCnt += 1
                    protCntb += 1
                else:
                    protCntb += 1
                    protein_dict[protName]['proteinIdentificationAs'] += [idAs]
                    protein_dict[protName]['proteinArea'] += areas
                    protein_dict[protName]['peptides'] += peptides
                    protein_dict[protName]['proteinCoverage'] += protCov
        else:
            pgCnt += 1
            for p in pg_dict[pg]['allMembers']:
                pgCnt2 += 1
                pn = p_dict[p]['desc']
                if "variant" in pn:
                    idAs = [pn.split(".")[1]]
                    protName = pn.split(".")[0]
                else:
                    protName = pn
                    idAs = ["ref"]
                others = pg_dict[pg]['allMembers']
                others.remove(p)
                otherNames = []
                for o in others:
                    otherNames.append(p_dict[o]['desc'])
                areas.append(p_dict[p]['area'])
                peptides.append(p_dict[p]['peptides'])
                protCov.append(p_dict[p]['coverage'])
                if protName in protein_dict:
                    prots.append(protName)
                    if 'proteinIdentification' not in protein_dict[protName]:
                        protein_dict[protName]['proteinIdentification'] = "putatively"
                        protein_dict[protName]['proteinIdentificationAs'] = [idAs]
                        protein_dict[protName]['proteinArea'] = areas
                        protein_dict[protName]['peptides'] = peptides
                        protein_dict[protName]['proteinCoverage'] = protCov
                        protein_dict[protName]['otherProteinsInGroup'] = otherNames
                        protCnt += 1
                        protCntb += 1
                    else:
                        protCntb += 1
                        protein_dict[protName]['proteinIdentificationAs'] += [idAs]
                        protein_dict[protName]['proteinArea'] += areas
                        protein_dict[protName]['peptides'] += peptides
                        protein_dict[protName]['proteinCoverage'] += protCov
                        if 'otherProteinsInGroup' not in protein_dict[protName]:
                            protein_dict[protName]['otherProteinsInGroup'] = otherNames
                        else:
                            protein_dict[protName]['otherProteinsInGroup'] += otherNames
    else:
        huoCnt += 1
print "processed protein groups representing more than one protein", pgCnt
print "processed protein groups representing the same protein", pgCntb
print "proteins in processed protein groups representing more than one protein", pgCnt2
print "different proteins matching db: ", len(set(prots))
print "protein groups without microbial proteins: ", huoCnt
print "retrieved proteomics results for ", protCntb, " proteins"
print "integrated proteomics results into ", protCnt, " gene entries"

# read file of essential genes and keep only best hits
header = 3
ess_dict = {}
essFile = "/work/projects/ecosystem_biology/MUST/CombinedAssembly/Annotation/hmms/" + FAM + "/" + LIB + "/contig.Prodigal.hmm.orfs.hits"
print "Reading essential gene information from ", essFile
ess_file = open(essFile, "r")
while 1:
    line = ess_file.readline()
    if line == "":
        break
    if header > 0:
        header -= 1
    else:
        line = line.rstrip()
        tab = line.split() # 0:target name, 1:accession, 2:query name, 3: accession, 4: E-value, 5: score, 6: bias, 7: E-value, 8:score, 9: bias, 10: exp, 11: reg, 12: clu, 13: ov, 14: env, 15: dom, 16: rep, 17: inc, 18: description of target
        if "#" not in tab[0]:
            if tab[0] not in ess_dict:
                ess_dict[tab[0]] = [tab[2],float(tab[4])]
            elif float(tab[4]) < ess_dict[tab[0]][1]:
                ess_dict[tab[0]] = [tab[2],float(tab[4])]
ess_file.close()

# update protein dictionary with essentiality information
for ess in ess_dict:
    if ess in protein_dict:
            protein_dict[ess]['essentialGene'] = ess_dict[ess][0] # I do not make an entry for unessential genes
ess_dict = {}

# read kegg annotation file for genes and update protein dictionary
header = 1
keggFile = "/work/projects/ecosystem_biology/MUST/CombinedAssembly/Annotation/hmms/" + FAM + "/" + LIB + "/contig.Prodigal.faakegg.besthits.tsv"
print "Reading KEGG annotations of genes from ", keggFile
kegg_file = open(keggFile, "r")
while 1:
    line = kegg_file.readline()
    if line == "":
        break
    if header == 1:
        header = 0
    else:
        line = line.rstrip()
        tab = line.split("\t") # 0:Gene, 1:KO, 2:maxScore, 3:hitNumber
        if tab[0] in protein_dict:
            kos = tab[1].split(";")
            protein_dict[tab[0]]['KO'] = kos # I do not make an entry for unannotated genes
kegg_file.close()

# read kegg network annotation file for genes and update protein dictionary
header = 1
keggNWFile = "/work/projects/ecosystem_biology/MUST/CombinedAssembly/Annotation/hmms/" + FAM + "/" + LIB + "/contig.Prodigal.faakegg.hitsNodes.tsv"
print "Reading KEGG network annotations of genes from ", keggNWFile
keggNW_file = open(keggNWFile, "r")
while 1:
    line = keggNW_file.readline()
    if line == "":
        break
    if header == 1:
        header = 0
    else:
        line = line.rstrip()
        tab = line.split("\t") # 0:Gene, 1:KO, 2:maxScore, 3:hitNumber
        if tab[0] in protein_dict:
            konws = tab[1].split(";")
            protein_dict[tab[0]]['node'] = konws # I do not make an entry for unannotated genes
keggNW_file.close()

# read metaCyc annotation file for genes and update protein dictionary
header = 1
metacycFile = "/work/projects/ecosystem_biology/MUST/CombinedAssembly/Annotation/hmms/" + FAM + "/" + LIB + "/contig.Prodigal.faametacyc.besthits.tsv"
print "Reading MetaCyc annotations of genes from ", metacycFile
metacyc_file = open(metacycFile, "r")
while 1:
    line = metacyc_file.readline()
    if line == "":
        break
    if header == 1:
        header = 0
    else:
        line = line.rstrip()
        tab = line.split("\t") # 0:Gene, 1:metaCycID, 2:maxScore, 3:hitNumber
        if tab[0] in protein_dict:
            mcs = tab[1].split(";")
            protein_dict[tab[0]]['metaCycID'] = mcs # I do not make an entry for unannotated genes
metacyc_file.close()

# read SwissProt annotation file for genes and update protein dictionary
header = 1
swissprotFile = "/work/projects/ecosystem_biology/MUST/CombinedAssembly/Annotation/hmms/" + FAM + "/" + LIB + "/contig.Prodigal.faaswissprot.besthits.tsv"
print "Reading Swiss-prot annotations of genes from ", swissprotFile
swissprot_file = open(swissprotFile, "r")
while 1:
    line = swissprot_file.readline()
    if line == "":
        break
    if header == 1:
        header = 0
    else:
        line = line.rstrip()
        tab = line.split("\t") # 0:Gene, 1:swissprotEC, 2:maxScore, 3:hitNumber
        if tab[0] in protein_dict:
            sps = tab[1].split(";")
            protein_dict[tab[0]]['swissprotEC'] = sps # I do not make an entry for unannotated genes
swissprot_file.close()

# read Pfam annotation file for genes and update protein dictionary
header = 1
pfamFile = "/work/projects/ecosystem_biology/MUST/CombinedAssembly/Annotation/hmms/" + FAM + "/" + LIB + "/contig.Prodigal.faapfam.besthits.tsv"
print "Reading Pfam annotations of genes from ", pfamFile
pfam_file = open(pfamFile, "r")
while 1:
    line = pfam_file.readline()
    if line == "":
        break
    if header == 1:
        header = 0
    else:
        line = line.rstrip()
        tab = line.split("\t") # 0:Gene, 1:pfamID, 2:maxScore, 3:hitNumber
        if tab[0] in protein_dict:
            pfs = tab[1].split(";")
            protein_dict[tab[0]]['pfamID'] = pfs # I do not make an entry for unannotated genes
pfam_file.close()

# read TIGR annotation file for genes and update protein dictionary
header = 1
tigrFile = "/work/projects/ecosystem_biology/MUST/CombinedAssembly/Annotation/hmms/" + FAM + "/" + LIB + "/contig.Prodigal.faatigrpfam.besthits.tsv"
print "Reading TIGR-Pfam annotations of genes from ", tigrFile
tigr_file = open(tigrFile, "r")
while 1:
    line = tigr_file.readline()
    if line == "":
        break
    if header == 1:
        header = 0
    else:
        line = line.rstrip()
        tab = line.split("\t") # 0:Gene, 1:tigrID, 2:maxScore, 3:hitNumber
        if tab[0] in protein_dict:
            tis = tab[1].split(";")
            protein_dict[tab[0]]['tigrID'] = tis # I do not make an entry for unannotated genes
tigr_file.close()

# read best annotation file for genes and update protein dictionary
header = 1
bestFile = "/work/projects/ecosystem_biology/MUST/CombinedAssembly/Annotation/hmms/" + FAM + "/" + LIB + "/besthitsAllDB.tsv"
print "Reading best hit functional annotations of genes from ", bestFile
best_file = open(bestFile, "r")
while 1:
    line = best_file.readline()
    if line == "":
        break
    if header == 1:
        header = 0
    else:
        line = line.rstrip()
        tab = line.split("\t") # 0:Gene, 1:ID, 2:maxScore, 3:hitNumber
        if tab[0] in protein_dict:
            best = tab[1].split(";")
            protein_dict[tab[0]]['bestAnnotation'] = best # I do not make an entry for unannotated genes
best_file.close()

# read Amphora2 annotation file for genes and update protein dictionary
header = 1
amphoraFile = "/work/projects/ecosystem_biology/MUST/CombinedAssembly/Annotation/hmms/" + FAM + "/" + LIB + "/amphora2.corrected.tsv"
print "Reading phylogenetic information of marker genes from ", amphoraFile
amphora_file = open(amphoraFile, "r")
while 1:
    line = amphora_file.readline()
    if line == "":
        break
    if header == 1:
        header = 0
    else:
        line = line.rstrip()
        tab = line.split("\t") # 0:gene, 1:marker, 2:kingdom, 3:phylum, 4:class, 5:order, 6:family, 7:genus, 8:species
        if tab[0] in protein_dict:
            king = []
            phy = []
            cla = []
            order = []
            fam = []
            gen = []
            spec = []
            protein_dict[tab[0]]['amphoraMarker'] = tab[1] #the names of the fields are different from R, because I cannot have "." in the keys
            king = tab[2].split("(")
            king[1] = king[1].replace(")","")
            protein_dict[tab[0]]['amphoraKingdom'] = [king[0],float(king[1])]
            if len(tab)>3:
                phy = tab[3].split("(")
                if len(phy)>1:
                    phy[-1] = phy[-1].replace(")","")
                    protein_dict[tab[0]]['amphoraPhylum'] = [phy[0],float(phy[-1])]
                else:
                    protein_dict[tab[0]]['amphoraPhylum'] = ["unassigned",0.0]
            if len(tab)>4:
                cla = tab[4].split("(")
                if len(cla)>1:
                    cla[-1] = cla[-1].replace(")","")
                    protein_dict[tab[0]]['amphoraClass'] = [cla[0],float(cla[-1])]
                else:
                    protein_dict[tab[0]]['amphoraClass'] = ["unassigned",0.0]
            if len(tab)>5:
                order = tab[5].split("(")
                if len(order)>1:
                    order[-1] = order[-1].replace(")","")
                    protein_dict[tab[0]]['amphoraOrder'] = [order[0],float(order[-1])]
                else:
                    protein_dict[tab[0]]['amphoraOrder'] = ["unassigned",0.0]
            if len(tab)>6:
                fam = tab[6].split("(")
                if len(fam)>1:
                    fam[-1] = fam[-1].replace(")","")
                    protein_dict[tab[0]]['amphoraFamily'] = [fam[0],float(fam[-1])]
                else:
                    protein_dict[tab[0]]['amphoraFamily'] = ["unassigned",0.0]
            if len(tab)>7:
                gen = tab[7].split("(")
                if len(gen)>1:
                    gen[-1] = gen[-1].replace(")","")
                    protein_dict[tab[0]]['amphoraGenus'] = [gen[0],float(gen[-1])]
                else:
                    protein_dict[tab[0]]['amphoraGenus'] = ["unassigned",0.0]
            if len(tab)>8:
                spec = tab[8].split("(")
                if len(spec)>1:
                    spec[-1] = spec[-1].replace(")","")
                    protein_dict[tab[0]]['amphoraSpecies'] = [spec[0],float(spec[-1])]
                else:
                    protein_dict[tab[0]]['amphoraSpecies'] = ["unassigned",0.0]
            # I do not make an entry for unannotated genes; the structure is different from how it used to be in R, as confidence levels are stored separately
amphora_file.close()

# the mOTU marker genes:
mOTUmarkers = ["COG0012","COG0016","COG0018","COG0172","COG0215","COG0495","COG0525","COG0533","COG0541","COG0552"]

# read best mOTU annotation files for genes and update protein dictionary
for mark in mOTUmarkers:
    header = 1
    motuBestFile = "/work/projects/ecosystem_biology/MUST/CombinedAssembly/Annotation/mOTUgenes/" + FAM + "/" + LIB + "/" + mark + ".bestHitPhylogeny.tsv"
    print "Reading phylogenetic information of marker gene", mark," from ", motuBestFile
    motuBest_file = open(motuBestFile, "r")
    while 1:
        line = motuBest_file.readline()
        if line == "":
            break
        if header == 1:
            header = 0
        else:
            line = line.rstrip()
            tab = line.split("\t")
            # 0:mOTU.species.annotation, 1:mOTU, 2:subjectID, 3:queryID, 4:percIdentity, 5:alignmentLength, 6:mismatches,
            ### 7:gapOpens, 8:q.start, 9:q.end, 10:s.start, 11:s.end, 12:evalue, 13:bitScore, 14:MarkerGene, 15:mOTU.linkage.Group.Name
            ### 16:Superkingdom, 17:Phylum, 18:Class, 19:Order, 20:Family, 21:Genus, 22:SpeciesCluster
            if tab[3] in protein_dict:
                protein_dict[tab[3]]['mOTUbestPercIdentity'] = float(tab[4]) #the names of the fields are different from R, because I cannot have "." in the keys
                protein_dict[tab[3]]['mOTUbestMarkerGene'] = tab[14]
                protein_dict[tab[3]]['mOTUbestSuperkingdom'] = tab[16]
                protein_dict[tab[3]]['mOTUbestPhylum'] = tab[17]
                protein_dict[tab[3]]['mOTUbestClass'] = tab[18]
                protein_dict[tab[3]]['mOTUbestOrder'] = tab[19]
                protein_dict[tab[3]]['mOTUbestFamily'] = tab[20]
                protein_dict[tab[3]]['mOTUbestGenus'] = tab[21]
                protein_dict[tab[3]]['mOTUbestSpeciesCluster'] = tab[22]
                protein_dict[tab[3]]['mOTUbestSpeciesAnnotation'] = tab[0] # I do not make an entry for unannotated genes;
    motuBest_file.close()
    
# read best present mOTU annotation files for genes and update protein dictionary
for mark in mOTUmarkers:
    header = 1
    motuBestFile = "/work/projects/ecosystem_biology/MUST/CombinedAssembly/Annotation/mOTUgenes/" + FAM + "/" + LIB + "/" + mark + ".bestPresentHitPhylogeny.tsv"
    print "Reading phylogenetic information of marker gene", mark,"of present mOTUs from ", motuBestFile
    motuBest_file = open(motuBestFile, "r")
    while 1:
        line = motuBest_file.readline()
        if line == "":
            break
        if header == 1:
            header = 0
        else:
            line = line.rstrip()
            tab = line.split("\t")
            # 0:mOTU.species.annotation, 1:mOTU, 2:subjectID, 3:queryID, 4:percIdentity, 5:alignmentLength, 6:mismatches,
            ### 7:gapOpens, 8:q.start, 9:q.end, 10:s.start, 11:s.end, 12:evalue, 13:bitScore, 14:MarkerGene, 15:mOTU.linkage.Group.Name
            ### 16:Superkingdom, 17:Phylum, 18:Class, 19:Order, 20:Family, 21:Genus, 22:SpeciesCluster
            if tab[3] in protein_dict:
                protein_dict[tab[3]]['mOTUpresentPercIdentity'] = float(tab[4]) # the names of the fields are different from R, because I cannot have "." in the keys
                protein_dict[tab[3]]['mOTUpresentMarkerGene'] = tab[14]
                protein_dict[tab[3]]['mOTUpresentSuperkingdom'] = tab[16]
                protein_dict[tab[3]]['mOTUpresentPhylum'] = tab[17]
                protein_dict[tab[3]]['mOTUpresentClass'] = tab[18]
                protein_dict[tab[3]]['mOTUpresentOrder'] = tab[19]
                protein_dict[tab[3]]['mOTUpresentFamily'] = tab[20]
                protein_dict[tab[3]]['mOTUpresentGenus'] = tab[21]
                protein_dict[tab[3]]['mOTUpresentSpeciesCluster'] = tab[22]
                protein_dict[tab[3]]['mOTUpresentSpeciesAnnotation'] = tab[0] # I do not make an entry for unannotated genes;
    motuBest_file.close()

#insert protein information into the contig dictionary
print "inserting information on ", len(protein_dict), "proteins into contig data"
contcount=0
for prot in protein_dict:
    protein_dict[prot]['gene'] = prot
    contig = re.sub(r'_[0-9]+$',"",prot)
    contig = contig.replace("gene","")
    if contig not in contig_dict:
        print "no contig information available for protein", prot, " on ", contig
    else:
        if "genes" not in contig_dict[contig]:
            contig_dict[contig]['genes'] = [protein_dict[prot]]
            contcount += 1
        else:
            contig_dict[contig]['genes'].append(protein_dict[prot])
protein_dict = {}
            
#insert RNA information into the contig dictionary
print "inserting information on ", len(rna_dict), "RNAs into contig data"
for rna in rna_dict:
    rna_dict[rna]['gene'] = rna
    contig = re.sub(r'_r[0-9]+$',"",rna)
    if contig not in contig_dict:
        print "no contig information available for RNA",rna, " on ", contig
    else:
        if "genes" not in contig_dict[contig]:
            contig_dict[contig]['genes'] = [rna_dict[rna]]
            contcount += 1
        else:
            contig_dict[contig]['genes'].append(rna_dict[rna])
rna_dict = {}

print "protein or RNA information for ", contcount, " contigs"

#make clusterInfo-like table
clusterStatFile = "/work/projects/ecosystem_biology/MUST/CombinedAssembly/Annotation/Bmaps/" + FAM + "/" + LIB + "/clusterInfoFromPy.tsv"
cluster_dict ={}
for cont in contig_dict:
    tmpCont = contig_dict[cont]
    complCnt = 0
    expCnt = 0
    expRead = 0
    proUCnt = 0
    proPCnt = 0
    ess = []
    ann = 0
    annKOL = 0
    annKONWL = 0
    annMCL = 0
    annSPL = 0
    annPfL = 0
    annTIL = 0
    rna16 = 0
    rna23 = 0
    annKO = []
    annKONW = []
    annMC = []
    annSP = []
    annPf = []
    annTI = []
    amPhy = []
    amGen = []
    mpPhy = []
    mpGen = []
    mp = []
    mbPhy = []
    mbGen = []
    mb = []
    geneno = 0
    if "genes" in tmpCont:
        tmpGenes = tmpCont['genes']
        geneno = len(tmpGenes)
        for g in tmpGenes:
            if g['completeness'] == "complete":
                complCnt += 1
            if g['aveCovRNAfw'] >= 0.01:
                expCnt += 1
            if g['readsRNAfw'] > 0:
                expRead += g['readsRNAfw']
            if 'proteinIdentification' in g:
                if g['proteinIdentification'] == "uniquely":
                   proUCnt += 1
                if g['proteinIdentification'] == "putatively":
                   proPCnt += 1
            if g['kind'] == "16S_rRNA":
                rna16 +=1
            if g['kind'] == "23S_rRNA":
                rna23 +=1
            if 'essentialGene' in g:
                ess.append(g['essentialGene'])
            if 'KO' in g or 'metaCycID' in g or 'swissProtEC' in g or 'pfamID' in g or 'tigrID' in g or 'bestAnnotation' in g:
                ann += 1
            if 'KO' in g:
                annKO += g['KO']
                annKOL += 1
            if 'node' in g:
                annKONW += g['node']
                annKONWL += 1
            if 'metaCycID' in g:
                annMC += g['metaCycID']
                annMCL += 1
            if 'swissprotEC' in g:
                annSP += g['swissprotEC']
                annSPL += 1
            if 'pfamID' in g:
                annPf += g['pfamID']
                annPfL += 1
            if 'tigrID' in g:
                annTI += g['tigrID']
                annTIL += 1
            if 'amphoraPhylum' in g:
                amPhy.append(g['amphoraPhylum'][0])
            if 'amphoraGenus' in g:
                amGen.append(g['amphoraGenus'][0])
            if 'mOTUbestMarkerGene' in g:
                mbPhy.append(g['mOTUbestPhylum'])
                mbGen.append(g['mOTUbestGenus'])
                mb.append(g['mOTUbestSpeciesAnnotation'])
            if 'mOTUpresentMarkerGene' in g:
                mpPhy.append(g['mOTUpresentPhylum'])
                mpGen.append(g['mOTUpresentGenus'])
                mp.append(g['mOTUpresentSpeciesAnnotation'])
    varlen = 0
    if "varPos" in tmpCont:
        varlen = len(tmpCont['varPos'])
    kP = []
    kG = []
    if "krakenPhylum" in tmpCont:
        kP = [tmpCont['krakenPhylum']]
    if "krakenGenus" in tmpCont:
        kG = [tmpCont['krakenGenus']]
    if tmpCont['cluster'] != "S":
        if tmpCont['cluster'] not in cluster_dict:
            cluster_dict[tmpCont['cluster']] = {'length' : tmpCont['length'],'contigs' : 1,'aveCov' : tmpCont['length']*tmpCont['aveCov'],
                                               'varPerMB' : varlen, 'genes' : geneno, 'RNA16' : rna16, 'RNA23' : rna23,
                                               'completeGenes' : complCnt, 'expressedGenes' : expCnt,'mappingRNAreads' : expRead, 'proteinsFoundUniquely': proUCnt,
                                               'proteinsFoundPutatively': proPCnt, 'ess' : ess,'annotated' : ann,
                                                'annotatedKO' : annKOL, 'annotatedKONW' : annKONWL, 'annotatedMetaCyc' : annMCL, 'annotatedSwissProt' : annSPL,
                                                'annotatedPfam' : annPfL, 'annotatedTIGR' : annTIL,
                                                'catsKO' : annKO, 'catsKONW' : annKONW, 'catsMetaCyc' : annMC, 'catsSwissProt' : annSP,
                                                'catsPfam' : annPf, 'catsTIGR' : annTI,
                                                'krakenPhyla' : kP,'krakenGenera' : kG,
                                                    'amphoraPhyla' : amPhy, 'amphoraGenera' : amGen,
                                                'motuPresentPhyla' : mpPhy, 'motuPresentGenera' : mpGen,
                                                'motuPresent' : mp, 'motuBestPhyla' : mbPhy, 'motuBestGenera' : mbGen, 'motuBest' : mb,
                                                'coordX' : [tmpCont['coords'][0]], 'coordY' : [tmpCont['coords'][1]]}
        else:
            cluster_dict[tmpCont['cluster']]['length'] += tmpCont['length']
            cluster_dict[tmpCont['cluster']]['contigs'] +=  1
            cluster_dict[tmpCont['cluster']]['aveCov'] += tmpCont['length']*tmpCont['aveCov']
            cluster_dict[tmpCont['cluster']]['varPerMB'] += varlen
            cluster_dict[tmpCont['cluster']]['genes'] += geneno
            cluster_dict[tmpCont['cluster']]['RNA16'] += rna16
            cluster_dict[tmpCont['cluster']]['RNA23'] += rna23
            cluster_dict[tmpCont['cluster']]['completeGenes'] += complCnt
            cluster_dict[tmpCont['cluster']]['expressedGenes'] += expCnt
            cluster_dict[tmpCont['cluster']]['mappingRNAreads'] += expRead
            cluster_dict[tmpCont['cluster']]['proteinsFoundUniquely'] += proUCnt
            cluster_dict[tmpCont['cluster']]['proteinsFoundPutatively'] += proPCnt
            cluster_dict[tmpCont['cluster']]['ess'] += ess
            cluster_dict[tmpCont['cluster']]['annotated'] += ann
            cluster_dict[tmpCont['cluster']]['annotatedKO'] += annKOL
            cluster_dict[tmpCont['cluster']]['annotatedKONW'] += annKONWL
            cluster_dict[tmpCont['cluster']]['annotatedMetaCyc'] += annMCL
            cluster_dict[tmpCont['cluster']]['annotatedSwissProt'] += annSPL
            cluster_dict[tmpCont['cluster']]['annotatedPfam'] += annPfL
            cluster_dict[tmpCont['cluster']]['annotatedTIGR'] += annTIL
            cluster_dict[tmpCont['cluster']]['catsKO'] += annKO
            cluster_dict[tmpCont['cluster']]['catsKONW'] += annKONW
            cluster_dict[tmpCont['cluster']]['catsMetaCyc'] += annMC
            cluster_dict[tmpCont['cluster']]['catsSwissProt'] += annSP
            cluster_dict[tmpCont['cluster']]['catsPfam'] += annPf
            cluster_dict[tmpCont['cluster']]['catsTIGR'] += annTI
            cluster_dict[tmpCont['cluster']]['krakenPhyla'] += kP
            cluster_dict[tmpCont['cluster']]['krakenGenera'] += kG
            cluster_dict[tmpCont['cluster']]['amphoraPhyla'] += amPhy
            cluster_dict[tmpCont['cluster']]['amphoraGenera'] += amGen
            cluster_dict[tmpCont['cluster']]['motuPresentPhyla'] += mpPhy
            cluster_dict[tmpCont['cluster']]['motuPresentGenera'] += mpGen
            cluster_dict[tmpCont['cluster']]['motuPresent'] += mp
            cluster_dict[tmpCont['cluster']]['motuBestPhyla'] += mbPhy
            cluster_dict[tmpCont['cluster']]['motuBestGenera'] += mbGen
            cluster_dict[tmpCont['cluster']]['motuBest'] += mb
            cluster_dict[tmpCont['cluster']]['coordX'].append(tmpCont['coords'][0]),
            cluster_dict[tmpCont['cluster']]['coordY'].append(tmpCont['coords'][1])
    else:
        if "S" not in cluster_dict:
            cluster_dict['S'] = {'length' : tmpCont['length'],'contigs' : 1,'aveCov' : tmpCont['length']*tmpCont['aveCov'],
                                  'varPerMB' : varlen, 'genes' : geneno, 'RNA16' : rna16, 'RNA23' : rna23,
                                  'completeGenes' : complCnt, 'expressedGenes' : expCnt,'mappingRNAreads' : expRead,'proteinsFoundUniquely': proUCnt,
                                  'proteinsFoundPutatively': proPCnt,
                                  'ess' : ess,'annotated' : ann,
                                  'annotatedKO' : annKOL, 'annotatedKONW' : annKONWL, 'annotatedMetaCyc' : annMCL, 'annotatedSwissProt' : annSPL,
                                  'annotatedPfam' : annPfL, 'annotatedTIGR' : annTIL,
                                  'catsKO' : annKO,'catsKONW' : annKONW, 'catsMetaCyc' : annMC, 'catsSwissProt' : annSP,
                                  'catsPfam' : annPf, 'catsTIGR' : annTI,
                                  'krakenPhyla' : kP,'krakenGenera' : kG,
                                  'amphoraPhyla' : amPhy, 'amphoraGenera' : amGen,
                                  'motuPresentPhyla' : mpPhy, 'motuPresentGenera' : mpGen,
                                  'motuPresent' : mp, 'motuBestPhyla' : mbPhy, 'motuBestGenera' : mbGen, 'motuBest' : mb,}
        else:
            cluster_dict['S']['length'] += tmpCont['length']
            cluster_dict['S']['contigs'] +=  1
            cluster_dict['S']['aveCov'] += tmpCont['length']*tmpCont['aveCov']
            cluster_dict['S']['varPerMB'] += varlen
            cluster_dict['S']['genes'] += geneno
            cluster_dict['S']['RNA16'] += rna16
            cluster_dict['S']['RNA23'] += rna23
            cluster_dict['S']['completeGenes'] += complCnt
            cluster_dict['S']['expressedGenes'] += expCnt
            cluster_dict['S']['mappingRNAreads'] += expRead
            cluster_dict['S']['proteinsFoundUniquely'] += proUCnt
            cluster_dict['S']['proteinsFoundPutatively'] += proPCnt
            cluster_dict['S']['ess'] += ess
            cluster_dict['S']['annotatedKO'] += annKOL
            cluster_dict['S']['annotatedKONW'] += annKONWL
            cluster_dict['S']['annotatedMetaCyc'] += annMCL
            cluster_dict['S']['annotatedSwissProt'] += annSPL
            cluster_dict['S']['annotatedPfam'] += annPfL
            cluster_dict['S']['annotatedTIGR'] += annTIL
            cluster_dict['S']['catsKO'] += annKO
            cluster_dict['S']['catsKONW'] += annKONW
            cluster_dict['S']['catsMetaCyc'] += annMC
            cluster_dict['S']['catsSwissProt'] += annSP
            cluster_dict['S']['catsPfam'] += annPf
            cluster_dict['S']['catsTIGR'] += annTI
            cluster_dict['S']['krakenPhyla'] += kP
            cluster_dict['S']['krakenGenera'] += kG
            cluster_dict['S']['amphoraPhyla'] += amPhy
            cluster_dict['S']['amphoraGenera'] += amGen
            cluster_dict['S']['motuPresentPhyla'] += mpPhy
            cluster_dict['S']['motuPresentGenera'] += mpGen
            cluster_dict['S']['motuPresent'] += mp
            cluster_dict['S']['motuBestPhyla'] += mbPhy
            cluster_dict['S']['motuBestGenera'] += mbGen
            cluster_dict['S']['motuBest'] += mb
print "collected data on ", len(cluster_dict), " clusters."
print "Writing cluster data to ", clusterStatFile
clusterStat_file = open(clusterStatFile, "w")
clusterStat_file.write("sample\tcluster\tlength\tcontigs\taveCov\tvarPerMB\tgenes\tcompleteGenes\texpressedGenes\tmappingRNAreads\tproteinsUniquely\tproteinsPutatively\t16SrRNAloci\t23SrRNAloci\tuniqueEss\tnumEss\t"+
                       "catsKO\tcatsMetaCyc\tcatsSwissProt\tcatsPfam\tcatsTIGR\tannotated\tannotatedKO\tannotatedMetaCyc\tannotatedSwissProt\tannotatedPfam\tannotatedTIGR\t" +
                       "phylaKraken\tgeneraKraken\tphylaAmphora\tgeneraAmphora\tphylaMotuPresent\tgeneraMotuPresent\tmotuPresent\tphylaMotuBest\tgeneraMotuBest\tmotuBest\tx\ty\n")
for clus in cluster_dict:
    length = str(cluster_dict[clus]['length'])
    contigs = str(cluster_dict[clus]['contigs'])
    aveCov = str(cluster_dict[clus]['aveCov']/cluster_dict[clus]['length'])
    varPerMB = str(cluster_dict[clus]['varPerMB']*1000000/cluster_dict[clus]['length'])
    genes = str(cluster_dict[clus]['genes'])
    complGenes = str(cluster_dict[clus]['completeGenes'])
    exprGenes = str(cluster_dict[clus]['expressedGenes'])
    exprReads = str(cluster_dict[clus]['mappingRNAreads'])
    uniqPros = str(cluster_dict[clus]['proteinsFoundUniquely'])
    putaPros = str(cluster_dict[clus]['proteinsFoundPutatively'])
    rnaSU = str(cluster_dict[clus]['RNA16'])
    rnaLU = str(cluster_dict[clus]['RNA23'])
    uniqueEss = str(len(set(cluster_dict[clus]['ess'])))
    numEss = str(len(cluster_dict[clus]['ess']))
    annotated = str(cluster_dict[clus]['annotated'])
    catsKO = str(len(set(cluster_dict[clus]['catsKO'])))
    annoKO = str(cluster_dict[clus]['annotatedKO'])
    catsKONW = str(len(set(cluster_dict[clus]['catsKONW'])))
    annoKONW = str(cluster_dict[clus]['annotatedKONW'])
    catsMetaCyc = str(len(set(cluster_dict[clus]['catsMetaCyc'])))
    annoMetaCyc = str(cluster_dict[clus]['annotatedMetaCyc'])
    catsSwissProt = str(len(set(cluster_dict[clus]['catsSwissProt'])))
    annoSwissProt = str(cluster_dict[clus]['annotatedSwissProt'])
    catsPfam = str(len(set(cluster_dict[clus]['catsPfam'])))
    annoPfam = str(cluster_dict[clus]['annotatedPfam'])
    catsTIGR = str(len(set(cluster_dict[clus]['catsTIGR'])))
    annoTIGR = str(cluster_dict[clus]['annotatedTIGR'])
    phylaKraken = enumerateTaxa(cluster_dict[clus]['krakenPhyla'])
    generaKraken = enumerateTaxa(cluster_dict[clus]['krakenGenera'])
    phylaAmphora = enumerateTaxa(cluster_dict[clus]['amphoraPhyla'])
    generaAmphora = enumerateTaxa(cluster_dict[clus]['amphoraGenera'])
    phylaMotuPresent = enumerateTaxa(cluster_dict[clus]['motuPresentPhyla'])
    generaMotuPresent = enumerateTaxa(cluster_dict[clus]['motuPresentGenera'])
    motuPresent = enumerateTaxa(cluster_dict[clus]['motuPresent'])
    phylaMotuBest = enumerateTaxa(cluster_dict[clus]['motuBestPhyla'])
    generaMotuBest = enumerateTaxa(cluster_dict[clus]['motuBestGenera'])
    motuBest = enumerateTaxa(cluster_dict[clus]['motuBest'])
    writeClusterList = [LIB,clus,length,contigs,aveCov,varPerMB,genes,complGenes,exprGenes,exprReads,uniqPros,putaPros,rnaSU,rnaLU,uniqueEss,numEss,
                        catsKO,catsKONW,catsMetaCyc,catsSwissProt,catsPfam,catsTIGR,annotated,annoKO,annoKONW,annoMetaCyc,annoSwissProt,annoPfam,annoTIGR,
                        phylaKraken,generaKraken,phylaAmphora,generaAmphora,phylaMotuPresent,generaMotuPresent,motuPresent,phylaMotuBest,generaMotuBest,motuBest,"NA\tNA"]
    if clus != "S":
        coordX = str(numpy.median(numpy.array(cluster_dict[clus]['coordX'])))
        coordY = str(numpy.median(numpy.array(cluster_dict[clus]['coordY'])))
        writeClusterList = [LIB,clus,length,contigs,aveCov,varPerMB,genes,complGenes,exprGenes,exprReads,uniqPros,putaPros,rnaSU,rnaLU,uniqueEss,numEss,
                            catsKO,catsKONW,catsMetaCyc,catsSwissProt,catsPfam,catsTIGR,annotated,annoKO,annoKONW,annoMetaCyc,annoSwissProt,annoPfam,annoTIGR,
                            phylaKraken,generaKraken,phylaAmphora,generaAmphora,phylaMotuPresent,generaMotuPresent,motuPresent,phylaMotuBest,generaMotuBest,motuBest,coordX,coordY]
    writeCluster = "\t".join(writeClusterList)
    clusterStat_file.write(writeCluster+"\n")
    cluster_dict[clus] = []
clusterStat_file.close()
cluster_dict = {}

#insert contig dictionaries into MongoDB
client = MongoClient()
db = client['mydb']
oldCollSize = db.must.count()
for cont in contig_dict:
    contig_dict[cont]['contig'] = cont
    db.must.insert_one(contig_dict[cont])
newCollSize = db.must.count()
print "contigs inserted into database:", newCollSize - oldCollSize
print "there are now", newCollSize, "documents in the collection."
