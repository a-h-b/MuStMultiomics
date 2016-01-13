In the MuSt, we used metagenomic/metatranscriptomic coverage for two main purposes: In the context of [binning contigs](automatic-clustering.md) into population-level genomes, we used average depth of metagenomic coverage of contigs. In the context of expression analysis, we used total numbers of fragments mapping to genes with functional annotations.

Both measures were inferred from mapping reads to genes or contigs. 

We mapped (filtered, trimmed) metagenomic reads to contigs using Bowtie2 (http://bowtie-bio.sourceforge.net/bowtie2/index.shtml).

```
bowtie2-build contigs.fa contigs
bowtie2 -x contigs -1 reads.screened.hg19.pair.1.fq -2 reads.screened.hg19.pair.2.fq -U reads.screened.hg19.single.fq -S DNAonContigs.sam -p 6 &>/DNAonContigs.log
```

Similarly, metagenomic reads can be mapped to genes.

```
bowtie2-build genes.rRNA.fa genes.rRNA
bowtie2 -x genes.rRNA -1 reads.screened.hg19.pair.1.fq -2 reads.screened.hg19.pair.2.fq -U reads.screened.hg19.single.fq -S DNAonGenesrRNA.sam -p 6 &>DNAonGenesrRNA.log
```

As we had stranded metatranscriptomic libraries, the mapping of RNA to reads looked a little different. To get reads mapping to the genes in sense (assuming, we use the index created above):

```
bowtie2 -x genes.rRNA --fr --nofw -1 reads.screened.hg19.pair.1.fq -2 reads.screened.hg19.pair.2.fq -U reads.screened.hg19.single.fq -S RNAonGenesrRNA.sam -p 6 &>RNAonGenesrRNA.log
```
To find all reads mapping antisense to genes could be done like this:

```
bowtie2 -x genes.rRNA --fr --norc -1 reads.screened.hg19.pair.1.fq -2 reads.screened.hg19.pair.2.fq -U reads.screened.hg19.single.fq -S revcompRNAonGenesrRNA.sam -p 6 &>revcompRNAonGenesrRNA.log
```

Now, to get the average depth of coverage (in this case of the contigs with metagenomic reads):

```
samtools view -bS DNAonContigs.sam > DNAonContigs.bam
samtools sort DNAonContigs.bam DNAonContigs.sorted
samtools depth DNAonContigs.sorted.bam > DNAonContigs.depth.txt
calculateCoverageAndGaps2.pl contigs.fa DNAonContigs.depth.txt > DNAonContigs.cov.tsv
```
The last line uses a custom script, which calculates and uses the length of the contig to calculate the average depth of coverage from the output of samtools' depth, because this omits nucleotides not covered at all.

To get the number of mapping (e.g. metatranscriptomic) fragments for each gene, we went to the Patrick May School of Wizadry and Oneliners (http://wwwen.uni.lu/lcsb/people/patrick_may)

```
grep -v -P "^\@" RNAonGenesrRNA.sam | cut -f 1,3 | sort | uniq | cut -f 2  | sort | uniq -c | perl -ane '$_=~/^\s+(\d+) (.+)$/;chomp($2); print "$2\t$1\n"; ' >RNAonGenesrRNA.mapped.tsv
```
This table can then be used to sum up the fragments mapping to genes with certain [functional annotations](functional-annotations.md) which can be used for differential analysis.

In addition to the described purposes, metagenomic or metatranscriptomic coverage of reads mapping to marker genes were used to infer community structures. For this purpose we applied the mOTU workflow (http://www.bork.embl.de/software/mOTU/) which forms part of MOCAT (http://vm-lux.embl.de/~kultima/MOCAT/index.html). Average depths of metagenomic or metatranscriptomic coverage were also used for general assessment of metagenomic and metatranscriptomic representations of genes, functions or contig bins within a sample. In very few instances, we compared these between samples after normalization to the total number of mapping nucleotides (sum of mapping lengths of all reads).


