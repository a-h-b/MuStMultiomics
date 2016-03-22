This workflow makes use of the functional annotation data and the completeness information stored in the [mongo database](mongo-database.md), collects the amino acid sequences of all complete phylogenetic marker gene predictions in all samples and builds a phylogenetic tree based on multiple sequence alignment for each class of marker genes. The information about membership of the marker genes in a contig [cluster](automatic-clustering.md) (binned reconstructed population-level genome) is retained, so closely related reconstructed genomes from different samples can be indentified as such.

In the first step, an R-script is run to retrieve all complete genes annotated as one of the mOTU [marker genes](annotate-phylogenetic-marker-genes.md) - COG0012, COG0016, COG0018, COG0172, COG0215, COG0495, COG0525, COG0533, COG0541, COG0552 - or __rpoB__ (TIGR02013) from the [mongo database](mongo-database.md). It returns a table with the gene IDs, the sample IDs and the contig cluster IDs of the complete marker genes for each sample and class of marker gene. (The script is documented [here](150816_getPhyloMarkers.R), but it depends on the specific file structure used in this project, so don't expect it to run anywhere else).

```
Rscript 150816_getPhyloMarkers.R
```

This table is then used as input to a perl script which extracts the genes of interest from the fasta file containing all gene predictions of a sample (names of samples are part of the combinedIds file). We append the output of the script to the same file, so we have one fasta file with all complete genes for each class of marker genes.

```
for lib in `cut -f 2 combinedIds`
do
  cd $lib/
  for marker in COG0012 COG0016 COG0018 COG0172 COG0215 COG0495 COG0525 COG0533 COG0541 COG0552
  do
    perl fastaProteinExtractAddSampleCluster.pl $marker.faa $marker.geneNamesClusters.tsv >> ../../trees/$marker.allCompleteGenes.faa
  done
  perl fastaProteinExtractAddSampleCluster.pl genes.faa rpoB.geneNamesClusters.tsv >> ../../trees/rpoB.allCompleteGenes.faa 
  cd ../
done
cd ../../trees/
```

Now, we use the [ETE toolkit](http://etetoolkit.org/) to build a phylogenetic tree for each of the marker gene classes.

```
for gene in rpoB COG0012 COG0016 COG0018 COG0172 COG0215 COG0495 COG0525 COG0541 COG0552 COG0533
do
  echo $gene
  python ete.py build --noimg --cpu 1 -a $gene.allCompleteGenes.faa -o allCompleteGenes.muscle_default-none-none-fasttree_default -w muscle_default-none-none-fasttree_default
done
```

The output of ETE includes a tree in Newick file format _.allCompleteGenes.faa.final_tree.nw_. We can read this into R using the package [ape](http://ape-package.ird.fr/). 
To find groups of marker genes which are close to each other, we can search the tree for polytomous points with more than 7 leaves.

```
library(ape)
library(geiger)

markers <- c("rpoB","COG0012","COG0016","COG0018","COG0172","COG0215","COG0495","COG0525","COG0533","COG0541","COG0552")
for(mark in markers){
  tall <- read.tree(paste(mark,".allCompleteGenes.faa.final_tree.nw",sep=""))
  degree <- tabulate(tall$edge[, 1])
  target <- which(degree > 7)
  
  for(group in target){
    write.table(data.frame("sample"=gsub("_.+","",tips(tall,group)),
                           "cluster"=gsub("M.+_","",gsub("_Contig.+","",tips(tall,group))),stringsAsFactors=F),
                           paste(mark,group,"clusters.tsv",sep="."),row.names=F,quote=F,sep="\t")
  }
}
```

This will write a tab-separated table for each group of closely related genes. The tables contain the ID of samples and the ID of the contig clusters (population-level genome) which contain a phylogenetic marker gene from the group of related genes. The files are named by a number that derives from the phylogenetic tree and the marker gene. 

Another useful package for phylogenetic tree visualization in R is [geiger](http://www.webpages.uidaho.edu/~lukeh/software.html). Using this and ape, we can visualize a tree for each marker gene and colour it according to the sample of origin. Or we can extract specific parts of the tree, like for example leaves formed by phylogenetic marker genes which belong to contig clusters with weak classification. Both examples are performed in the script [`150819_MUST_tree.R`](150819_MUST_tree.R), which contains some nice functions for visualization. The script is however based on the MuSt sample IDs and MuSt file structure and should not be run as is on other data sets.
