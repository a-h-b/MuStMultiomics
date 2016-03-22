To annotate the predictions of protein coding genes with functions, we used different sets of HMM databases (for KEGG KOs, enzymes from MetaCyc and Uniprot, PFAM families, and TIGRFAM families). The databases for Pfam and TIGRFAM were downloaded from the [Pfam](ftp://ftp.ebi.ac.uk/pub/databases/Pfam/) and [TIGRFAM](ftp://ftp.jcvi.org/pub/data/TIGRFAMs/) websites. The other databases were in-house, as decribed [here](http://pubs.rsc.org/en/Content/ArticleLanding/2009/MB/b915913b#!divAbstract). 

HMMER 3.1 (http://hmmer.janelia.org/) was run on each database to find hits ( _contig.Prodigal.faa_ is the file with amino acid sequences of the protein coding gene predictions).

```
hmmsearch --cpu 12 --noali --notextw --tblout contig.Prodigal.faa.kegg.hmmscan KO.hmm contig.Prodigal.faa >/dev/null
hmmsearch --cpu 12 --noali --notextw --tblout contig.Prodigal.faa.metacyc.hmmscan metacyc.hmm contig.Prodigal.faa >/dev/null
hmmsearch --cpu 12 --noali --notextw --tblout contig.Prodigal.faa.Pfam-A.hmmscan Pfam-A.hmm contig.Prodigal.faa >/dev/null
hmmsearch --cpu 12 --noali --notextw --tblout contig.Prodigal.faa.tigrpfam.hmmscan tigrpfam.hmm contig.Prodigal.faa >/dev/null
hmmsearch --cpu 12 --noali --notextw --tblout contig.Prodigal.faa.swissprot.hmmscan swissprot.hmm contig.Prodigal.faa >/dev/null
```
The output from HMMER was parsed using a [script](consolidate_hmmscan_results.pl) which returns all hits in a better format (tab-separated).

```
perl consolidate_hmmscan_results.pl contig.Prodigal.faa contig.Prodigal.faa.kegg.hmmscan contig.Prodigal.faa.metacyc.hmmscan contig.Prodigal.faa.Pfam-A.hmmscan contig.Prodigal.faa.swissprot.hmmscan contig.Prodigal.faa.tigrpfam.hmmscan
```
There is a [second version for this script](consolidate_hmmscan_results_justKEGG.pl) which does exactly the same, but only for the KO annotations.
The formatted output is then parsed once more to retain only the best annotation out of all databases, or to get the top or top x percent of best hits from any database. To get the best annotation per gene, we can use [`150310_MUST_hmmBestAll.py`](150310_MUST_hmmBestAll.py).

```
python 150310_MUST_hmmBestAll.py contig.Prodigal.faakegg.tsv contig.Prodigal.faametacyc.tsv contig.Prodigal.faaswissprot.tsv contig.Prodigal.faapfam.tsv contig.Prodigal.faatigrpfam.tsv -g $(grep ">" contig.Prodigal.faa | wc -l)
```
Alternatively, or additionally (in the case of MuSt), we can get [the best annotation per database](150705_MUST_hmmParse.py):

```
python 150705_MUST_hmmParsePfam.py contig.Prodigal.faapfam.tsv pfamID -g $(grep ">" contig.Prodigal.faa | wc -l) -k
python 150705_MUST_hmmParse.py contig.Prodigal.faametacyc.tsv metaCycID -g $(grep ">" contig.Prodigal.faa | wc -l)
python 150705_MUST_hmmParse.py contig.Prodigal.faaswissprot.tsv swissprotEC -g $(grep ">" contig.Prodigal.faa | wc -l)
python 150705_MUST_hmmParse.py contig.Prodigal.faatigrpfam.tsv tigrID -g $(grep ">" contig.Prodigal.faa | wc -l)
python 150705_MUST_hmmParse.py contig.Prodigal.faakegg.tsv KO -g $(grep ">" contig.Prodigal.faa | wc -l)
```
The IDs in the Pfam database have a slightly different format, so there is an [extra script](150705_MUST_hmmParsePfam.py) for this.
As a last alternative, we can annotate genes with the nodes of a [reconstructed metabolic network](reconstructed-KO-network) based on the KO annotations. The [script for this](150705_MUST_keggParseNW.py) uses a table (provided here in [`150705_KOs_in_NW.tsv`](150705_KOs_in_NW.tsv)) with the link between nodes and KOs (which node contains which KOs) which is built by the network reconstruction script.

```
python 150705_MUST_keggParseNW.py contig.Prodigal.faakegg.tsv -g $(grep ">" contig.Prodigal.faa | wc -l)
```
As you can see, I always use the -g option of these scripts, which sets a sample-specific threshold for annotations to be accepted. In accordance with the HMMer manual, this is log2 of the number of predicted genes.

Now we have the annotations for each gene, but we might also like to know how many reads (representing how many fragments) map to the genes with a functional annotation to do some expression analysis. This can be done using very similar a family of [scripts](150322_bestHmmReadParse.py), which also take a gene-wise read [coverage](calculating-coverage.md) table as input. In the MuSt, the gene expression analysis was based on the best annotation out of all databases:

```
python 150322_bestHmmReadParse.py contig.Prodigal.faakegg.tsv contig.Prodigal.faametacyc.tsv contig.Prodigal.faaswissprot.tsv contig.Prodigal.faapfam.tsv contig.Prodigal.faatigrpfam.tsv DNAonGenesrRNA.cov.tsv -g $(grep ">" contig.Prodigal.faa | wc -l)
```
There is also a version of this script which calculates the average coverage depth of each function from all genes annotated with this function ([`150415_bestHmmAveCovParse.py`](150415_bestHmmAveCovParse.py)). The results of this were not used within the MuSt.

```
#python 150415_bestHmmAveCovParse.py contig.Prodigal.faakegg.tsv contig.Prodigal.faametacyc.tsv contig.Prodigal.faaswissprot.tsv contig.Prodigal.faapfam.tsv contig.Prodigal.faatigrpfam.tsv DNAonGenesrRNA.cov.tsv -g $(grep ">" contig.Prodigal.faa | wc -l)
```
Alternatively, only single databases can be used. In the case of MuSt, we [used](150630_keggReadParse.py) the KO database.

```
python 150630_keggReadParse.py contig.Prodigal.faakegg.tsv DNAonGenesrRNA.mapped.tsv -g $(grep ">" contig.Prodigal.faa | wc -l) -p 100.0
```


