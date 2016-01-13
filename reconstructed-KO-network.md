In the MuSt, we used a previously established method (http://www.nature.com/articles/npjbiofilms20157) to reconstruct a community-wide metabolic network from KEGG orthologous groups (KOs) by forming edges between KO nodes, if the KOs share reactants. The reconstruction essentially uses the same script as in the previous project (INSERT LINK), with only one significant change in the combination of KOs in nodes. (In both projects, KOs of enzymes which have exactly the same reactants are combined into one nodes. The rationale behind this is that KOs with the exact same reactants are often subunits of the same enzyme. Within the concept of the network reconstructions and the topological analyses, these cases are better represented by single nodes than tight clusters of KOs.) In the MuSt workflow, the KOs are only combined into nodes, if they have adjacent KO numbers. This is due to the fact that KOs forming the subunits of one enzyme usually have adjacent numbers and KOs with the exact same reactants but non-adjacent numbers are often iso-enzymes or even enzymes catalyzing inverse reactions. To use the script `140630_MUST_NW.R`, several tables have to be downloaded from http://rest.kegg.jp:
http://rest.kegg.jp/link/Ko/pathway
http://rest.kegg.jp/link/Ko/rn
http://rest.kegg.jp/link/rn/rp
http://rest.kegg.jp/link/rp/rc/
http://rest.kegg.jp/link/rc/rn
http://rest.kegg.jp/list/rp
http://rest.kegg.jp/list/cpd
We also use the file ko2des_clean.txt which contains the descriptions for every KO in the HMM database.

The script gives out the reconstructed network, as well as a file with all KOs and nodes in the network. This can be used to find the expression levels of genes [functionally annotated](functional-annotations.md) with these KOs and restrict the network to only nodes represented by (expressed) genes.

Based on expression levels, sub-modules with significant differences under different conditions (such as samples from different individuals, families, dependent on disease etc) can be described. We used the R-package BioNet (https://bioconductor.org/packages/release/bioc/html/BioNet.html) for this. BioNet can be applied using the program Heinz (https://github.com/ls-cwi/heinz) or an R-implementation (FastHeinz). Function to either use the pure-R version of BioNet with output from DESeq2 (https://bioconductor.org/packages/release/bioc/html/DESeq2.html) or call Heinz on the command line from R (using `runHeinz.sh`) can be found in `140630_MUST_NW.R`.
These functions supplied with BioNet to make a graphical representation of the highest scoring module.

To make a more fancy plot, which also displays the relative gene/transcript/protein copy numbers linked to the binned population-level genome reconstructions, further functions can be found in `plotModules_omicLevels.R`. These functions, in turn, use data retrieved from a [mongo database](mongo-database.md). These plots are very rich in specific information, so the script should not be expected to work without this information.


