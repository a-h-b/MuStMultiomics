For the MuSt Multiomics analysis, we built sample specific search databases. To do this we used the genome sequences of the probands who donated the faecal samples to personalize the human proteins in the search databases. On the other hand, we used the genes predicted from the assemblies of metagenomic and metatranscriptomic reads of each sample. However, assemblies don't reflect the strain level variation encountered in a microbial community but rather a consensus. Therefore we also called variants. First, metagenomic and metatranscriptomic reads were mapped to the contigs (see [calculating coverage](calculating-coverage.md)). The names of the resulting .bam-files were saved in `bams.DNARNA.txt`. We then used Platypus (http://www.well.ox.ac.uk/platypus; version 0.7.9.1) to call variants:

```
Platypus.py callVariants --refFile=contigs.fa --bamFiles=bams.DNARNA.txt --nCPU=12 -o DNARNAonContigs.vcf
```
bgzip and tabix (https://github.com/samtools/htslib) were used to compress the .vcf-files.

```
bgzip DNARNAonContigs.vcf
tabix DNARNAonContigs.vcf.gz
```
vcftools (https://vcftools.github.io/man_latest.html) was used to filter the variants:

```
vcftools --gzvcf DNARNAonContigs.vcf.gz --remove-filtered GOF --remove-filtered QD --remove-filtered hp10 --remove-filtered Q20 --remove-filtered MQ --remove-filtered SC --remove-filtered QualDepth --remove-filtered badReads --recode --stdout | bgzip > DNARNAonContigsFiltered.vcf.gz
bgzip -d DNARNAonContigsFiltered.vcf.gz
```

For the next step, we changed some code written by Nic Pinel (http://www.eafit.edu.co/docentes-investigadores/Paginas/nicolas-pinel.aspx) for his project with Emilie Muller (http://www.nature.com/ncomms/2014/141126/ncomms6603/full/ncomms6603.html). Beside the .vcf this [script](variant_annotateRepairedTab.pl) uses the .tab output of prodigal (https://github.com/hyattpd/prodigal/releases/) which had been used to call the genes and the assembled contigs (which were formatted into blocks of width 80 using the fasta_formatter of the fastx suit (http://hannonlab.cshl.edu/fastx_toolkit/).

```
perl variant_annotateRepairedTab.pl -v DNARNAonContigsFiltered.vcf -a gene.tab -s contigs.formatted.fa -p
```
This  [script](variant_annotateRepairedTab.pl) gives out two fasta files with the nucleotide and amino acid sequences of all predicted proteins, including the original predictions and the variants, as well as a .tab file with the variant gene positions in the prodigal-style. The resulting .tab-file was concatenated with the original tab file. This was used to [remove](trypsinStartEnd.pl) incomplete protein fragments which would not form tryptic peptides.

```
cat gene.tab contigs.formatted.variants.tab >> genes.variants.tab
perl trypsinStartEnd.pl contigs.formatted.variants.faa genes.variants.tab >> genes.variants.endsRemoved.faa
```

Before the microbial protein predictions were concatenated with the human proteins, they were [renamed](rename4proteomics.pl) to have a unique ID starting with sp|1xxxxxxx|.

```
perl rename4proteomics.pl genes.variants.endsRemoved.faa 1 >> genes.variants.endsRemoved.Renamed.faa
```

