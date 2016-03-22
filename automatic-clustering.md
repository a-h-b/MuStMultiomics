One of the challenges in the MuSt was to link genes with functions of interest to their genomic contexts. One option to achieve this is reference-independent binning. During his time in the Eco-Systems Biology group at LCSB, [Cédric Laczny](https://github.com/claczny?tab=repositories) developed [VizBin](http://claczny.github.io/VizBin/), which uses Barnes-Hut Stochastic Neighbour Embedding of high-dimensional _kmer_ frequency data to build 2D maps of contigs. Within VizBin, clusters within these maps can be nicely picked manually. However, in the MuSt we wanted to pick as many clusters as possible from each sample and we had 36 samples with multi-omic data (and 53 with metagenomic assemblies), so manual picking really was no option. In addition, we wanted to have the cluster picking informed by metagenomic coverage and the presence of single-copy essential genes, which was not possible with VizBin back then (it is now!). 
The essential single-copy genes, as compiled by [Dupont et al, 2011](http://www.nature.com/ismej/journal/v6/n6/full/ismej2011189a.html) were called using HMMs (https://github.com/MadsAlbertsen/multi-metagenome/blob/master/R.data.generation/essential.hmm) provided by [Mads Albertsen](http://madsalbertsen.github.io/multi-metagenome/) using a wrapper script for HMMER (http://hmmer.janelia.org/).

```
hmmsearch --tblout hmm.orfs.hits --cut_tc --notextw essential.hmm contig.Prodigal.faa > hmmer.out
```
(If this does not work for you, it may be because you have spaces in the fasta headers of your gene predictions. You can do `cut -f1 -d" " contig.Prodigal.faa > contig.Prodigal.renamed.faa` to help yourself before running HMMER on the renamed-file).

There is one further trick we played here: Cédric had noticed quite early on in working with VizBin and [its predecessors](http://www.nature.com/articles/srep04516) that 16S and 23S rRNAs are so strongly conserved that they do not cluster with their genomes, unless the flanking regions are long enough to significantly impact on the genomic signature of the contig. We therefore removed 16S and 23S genes as much as possible without cutting contigs to lengths below 1000nt (which we need for a good signature). This cutting is done using .gff outputs from Barrnap (http://www.vicbioinformatics.com/software.barrnap.shtml) and a [perl script](fastaExtractCutRibosomal1000.pl) (using Bioperl (http://bioperl.org/), so we format the contigs first using fastx's fasta_formatter (http://hannonlab.cshl.edu/fastx_toolkit/) to avoid lines longer than 2^15 characters). The resulting fasta file includes all contigs with at least 1000 nt length. Names are conserved.

```
fasta_formatter -i contigs.fa -o contigs.formatted.fa -w 80
fastaExtractCutRibosomal1000.pl contigs.formatted.fa rRNAgenes.euk.gff rRNAgenes.arc.gff rRNAgenes.bac.gff rRNAgenes.mito.gff contigs.1000.rRNAcut.fa rRNAcutting.log
```
Then, we ran the _kmer_ frequency analysis on the cut contigs and the dimension reduction as implemented in VizBin non-interactively. As the coordinates of the 2D maps are returned in order, but without a name, a file with the contig names is also generated.

```
bash runBHSNE.sh contigs.1000.rRNAcut.fa .
grep ">" contigs.1000.rRNAcut.fa > contigs.1000.rRNAcut.names.txt
```

We then used DBSCAN (http://www.dbs.ifi.lmu.de/Publikationen/Papers/KDD-96.final.frame.pdf), as implemented in R (https://cran.r-project.org/web/packages/fpc/index.html) to bin clusters of contigs. Afterwards, we analyze the number and uniqueness of single-copy essential genes within the bins. If essential genes are present (more or less) in single copies, the clusters are accepted. If there are duplicates of the essential genes, we use the [average metagenomic coverage depth](calculating-coverage.md) of the contigs to separate the contigs into sub-bins with unimodal coverage distributions. The [script](autoCluster.R) uses Hartigans' diptest as implemented in the R-package diptest (https://cran.r-project.org/web/packages/diptest/index.html) to test for unimodality and a normal mixture model from the R-package mixtools (https://cran.r-project.org/web/packages/mixtools/index.html) to find a cut-off between groups of contigs with unimodal, (on a log scale) normally distributed coverage. Sub-bins are re-analyzed for uniqueness of the single-copy essential genes and the whole process is repeated three times with increasing stringency applied to the reachable points. All of these steps are performed by the script [`autoCluster.R`](autoCluster.R).

```
Rscript autoCluster.R 
```
This script uses an R workspace which we created previously and which contains all the information used during the binning. [`makeWSvarAnnoCorrect.R`](makeWSvarAnnoCorrect.R) documents how this workspace was built. I don't recommend running this script, as it is full of paths to files which will only be found on our system. 

The outputs of `autoCluster.R` include several files which document the cut-offs and reachability estimates used in the different iterations. Most importantly, a table with the cluster/bin membership of every contig is produced ( _contigs2clusters.tsv_ ). In addition, a 2D BH-SNE map with the final bins is created and for the relatively well recovered genomes (bins with >67% of the single-copy essential genes) bed files for the contigs and genes are created, as well as a visual evaluation of the length, coverage, taxonomy and functional content of the genomes.

The autoCluster method is being developed further. A version to be used with the output of the integrated metagenomic and metatranscriptomic assembly and analysis pipeline IMP, developed by Shaman Narayanasamy and Yohan Jarosz is available.


