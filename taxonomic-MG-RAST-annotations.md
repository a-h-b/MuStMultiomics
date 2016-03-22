MG-RAST (http://metagenomics.anl.gov/) can be used in many ways to analyze and annotate metagenomic (or metatranscriptomic) data. Independent of the kind of sequences that were uploaded, annotations are achieved using sequence alignment via Blast. Annotations can be accessed using the MG-RAST API (http://api.metagenomics.anl.gov/api.html) or the "Download" page.
Within the MuSt project, we analyzed both functional and taxonomic annotations. However, we found that functional annotation can be improved by [using HMMs](functional-annotations). In addition, we used the genomic context of the contigs and binning of the contigs to establish taxonomy of prokaryotic (mostly bacterial) genes. Therefore, the taxonomic assignments of potential eukaryotic and viral genes were the main application of MG-RAST within the MuSt.
So we uploaded the gene predictions (protein coding and bacterial rRNAs). If you want to check out our data, you can find the accession numbers in MGRASTaccessions.tsv. We then downloaded the RefSeq-based organism annotations via the Download-page. The downloaded tables were mined for eukaryotic genes (except for remaining human sequences) and viral genes (except the phiX174 genome, which is spiked into illumina sequencing libraries) using the R-script `MGRASTgeneLevelTax.R`. As one gene can be aligned to genes of several taxa, we kept the assignment with the highest bitscore for each gene. In addition, we kept all assignments with a score that was at least 80 % of the highest bitscore for this gene. To find the lowest common ancestor (LCA) and the species, genus, family, order, class, phylum and kingdom level taxonomy of each gene, the script calls the script `parse_taxbrowser_MGRAST.py`. This script relies heavily on a parser written by Romain Studer (http://evosite3d.blogspot.de/2013/06/browsing-ncbi-taxonomy-with-python.html) and the NCBI taxonomy. The output of the R-script is a table for each sample with the eukaryotic genes and a table with the viral genes (and some temporary files). Each table contains the gene name, the LCA, and the names at species, genus, family, order, class, phylum and kingdom level.
All additional information about these genes, such as coverage or detection at the protein level, were retrieved from the [Mongo database](mongo-database.md).
