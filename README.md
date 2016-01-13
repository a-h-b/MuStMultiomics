This repository contains code used in the multiomic analyses of faecal microbiota from four families with several cases of T1DM ( __MuSt__ ).

to build a search data base for [proteomics](proteomics-data-base) from predicted proteins and their variants:
rename4proteomics.pl
trypsinStartEnd.pl
variant_annotateRepairedTab.pl
variants_annotateTab4Stats.pl
variants_locateType.pl

to parse [functional annotations](functional-annotations) of gene predictions (some including [coverage](calculating-coverage)):
150310_MUST_hmmBestAll.py
150705_MUST_hmmParse.py
150705_MUST_hmmParsePfam.py
consolidate_hmmscan_results.pl
consolidate_hmmscan_results_justKEGG.pl
150705_MUST_keggParseNW.py
ko2des_clean.txt
calculateCoverageAndGaps2.pl
150322_bestHmmReadParse.py
150415_bestHmmAveCovParse.py
150630_keggReadParse.py

to annotate phylogenetic [marker genes with the taxonomy](annotate-phylogenetic-marker-genes) of the best hit from the mOTU database:
fastaExtractWithCoordBase1.pl
getValidGenes1.R
getValidGenes2new.R
getHitPhylogenyNew.R

to parse taxonomy of [MG-RAST annotations](taxonomic-MG-RAST-annotations) of genes:
parse_taxbrowser_MGRAST.py
MGRASTgeneLevelTax.R

to automatically [cluster](automatic-clustering) contigs based on nucleotide signature (BH-SNE maps), DNA coverage and essential genes:
autoCluster.R
fastaExtractCutRibosomal1000.pl
makeWSvarAnnoCorrect.R

to [reconstruct](reconstructed-KO-network) a metabolic network from KOs and analyse it:
140630_MUST_NW.R
with file 150705_KOs_in_NW.tsv
runHeinz.sh
plotModules_omicLevels.R

to feed a [mongo database](mongo-database) with all the data from MuSt and retrieve some of the data:
150928_mongifyMust.py
151020_funOIMongoWS.R
150928_MUST_relatedClusterWSFromMongo.R
eukaryoticGenesMongo.R
virusGenesMongo.R
getNWexprMongoAllSamples.R