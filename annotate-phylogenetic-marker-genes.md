For every assembly, phylogenetic marker genes were called from the gene predictions using fetchMG.pl (http://vm-lux.embl.de/~mende/fetchMG/cite.html).

```
perl fetchMG.pl -o marker_genes -t 4 -d contig.Prodigal.fna  -m extraction contig.Prodigal.faa
```
The mOTU software comes with a collection of phylogenetic marker genes padded with their genomic context ( _mOTU.v1.padded_ ). To annotate the phylogenetic marker genes found in our assembly, we aligned our genes with the genes in this collection belonging to mOTUs that had been detected in the sample of interest in the read level analysis (standard mOTU workflow) using blast, took the taxonomic annotation of the best hit and gave it to our gene. We call this the best present hit. In a second, alternative workflow, we took the taxonomic annotation of the best hit overall, without the limitation to the previously detected mOTUs. Here is the detailed description of how this was done.

Beside the marker genes collection, the mOTU software uses a mapping file _mOTU.v1.padded.motu.map_ which contains the gene Sequence ID, the length???, the COG ID and an internal mOTU ID. There is also the mapping file _mOTU.v1.map.txt_ containing the Sequence ID, the COG ID, the internal marker mOTU ID, and the mOTU linkage Group Name. In addition, _mOTU.v1.padded.motu.linkage.map_ contains the internal marker mOTU ID, and the taxonomy on the following levels: Domain, Phylum, Class, Order, Family, Genus, and Species, and, if applicable, the Species cluster, the Species cluster annotation and/or the mOTU species annotation (only the last column contains both Cluster and motu_linkage_group annotations).

We set up a directory for each COG present in the mapping files. We then put the Sequence IDs from the mapping file _mOTU.v1.map.txt_ into the directory of the COG they belong to.

```
for cog in `cut -f 3 mOTU.v1.padded.motu.map | sort | uniq`; do mkdir $cog;done
for cog in `cut -f 3 mOTU.v1.padded.motu.map | sort | uniq`; do grep $cog mOTU.v1.map.txt | cut -f 1 > $cog/$cog.mOTU.v1.map.txt ;done
```
There is also a file with the coordinates of the marker genes within the padded gene collection called _mOTU.v1.padded.coord_. From this extracted the coordinates for each gene and placed it in a file in the respective COG directory. The coordinates were then used to extract the gene sequences (without the padding) from the main database and place a fasta with all sequences belonging to one COG in the respective directory.

```
for cog in `cut -f 3 mOTU.v1.padded.motu.map | sort | uniq`; do for gene in $(cat $cog/$cog.mOTU.v1.map.txt); do grep $gene -w mOTU.v1.padded.coord >> $cog/$cog.mOTU.v1.padded.coord ;done; done
for cog in `cut -f 3 mOTU.v1.padded.motu.map | sort | uniq`; do cd $cog/; perl fastaExtractWithCoordBase1.pl ../mOTU.v1.padded $cog.mOTU.v1.padded.coord >> $cog.mOTU.v1.padded; cd ../; done
```

We also set up a directory for every sample, into which we copied the nucleotide sequences of the phylogenetic marker genes called on our assemblies (one file per COG).

```
for cog in `cut -f 3 mOTU.v1.padded.motu.map | sort | uniq`; do cp marker_genes/$cog.fna . ; done
```
This could already be used to find the best hit out of all mOTU marker genes, but we wanted to make use of the read-level analysis and limit the search only to those mOTUs that had been detected there. Therefore, we went into the results of the mOTU analysis, took the table _mOTU.counts.gz_, unpacked it, and extracted the mOTUs that were found. As the mOTU analysis had been done for metagenomic and metatranscriptomic reads, we combined both lists.

```
ln -s RESULTS/mOTU.counts.gz metaG.mOTU.counts.gz
gzip -d -c metaG.mOTU.counts.gz > metaG.mOTU.counts	
grep motus.processing -v metaG.mOTU.counts | cut -f 1 > mOTUs.present
ln -s RESULTS/mOTU.counts.gz metaT.mOTU.counts.gz
gzip -d -c metaT.mOTU.counts.gz > metaT.mOTU.counts
grep motus.processing -v metaT.mOTU.counts | cut -f 1 >> mOTUs.present
cat mOTUs.present | sort | uniq > mOTUs.present.all
```

With the help of the R script [`getValidGenes1.R`](getValidGenes1.), we then searched the mapping files _mOTU.v1.padded.motu.linkage.map_ and _mOTU.v1.map.txt_ for the genes belonging to the mOTUs found in the sample. The genes were then extracted from the marker gene collection in COG-wise fashion into a new fasta file.

```
for cog in `cut -f 3 ../../mOTU.v1.padded.motu.map | sort | uniq` 
do
	Rscript ../../getValidGenes1.R $cog
	perl testFastaExtract.pl ../../$cog/$cog.mOTU.v1.padded $cog.present.mOTU.genenames.txt > $cog.present.genes.fna
done
```
>(The leading ../../ in this example come from the directory structure with a directory per sample (ordered by family). They are here because they are also in the R script.

We then make a separate blast database with the limited genes for each COG and perform blastn for every gene annotated as this COG in our assembly. Already here, we limit the result to the best hit.

```
for cog in `cut -f 3 mOTU.v1.padded.motu.map | sort | uniq` 
do
	makeblastdb -in $cog.present.genes.fna -dbtype nucl -parse_seqids
	blastn -db $cog.present.genes.fna -query $cog.fna -outfmt 7 -out $cog.bestPresentHits.tsv -max_target_seqs 1
done
```
The annotations of the best hit on the different taxonomic levels are then taken from _mOTU-LG.v1.annotations.txt_, _mOTU.v1.padded.motu.linkage.map_ and _mOTU.v1.map.txt_ using the R script `getValidGenes2new.R`. 


In the alternative workflow, which finds the best hit out of the whole marker gene collection, we just built a blast database with all marker genes (per COG). 

```
for cog in `cut -f 3 mOTU.v1.padded.motu.map | sort | uniq`
do 
	cd $cog/
	makeblastdb -in $cog.mOTU.v1.padded -dbtype nucl -parse_seqids 
	cd ../
done
```
We then performed a blastn search for each sample using the predicted genes with this COG annotation from our assemblies. The taxonomy of this hit was again filtered from _mOTU-LG.v1.annotations.txt_, _mOTU.v1.padded.motu.linkage.map_ and _mOTU.v1.map.txt_ , this time using the script [`getHitPhylogenyNew.R`](getHitPhylogenyNew.R).

```
for cog in `cut -f 3 ../../mOTU.v1.padded.motu.map | sort | uniq`
do
	blastn -db ../../$cog/$cog.mOTU.v1.padded -query $cog.fna -outfmt 7 -out $cog.bestHits.tsv -max_target_seqs 1
	Rscript ../../getHitPhylogenyNew.R $cog
done
```

Finally, we gathered all annotations of marker genes within the binned population-level genomes. We found that in most cases, annotation was unanimous.


