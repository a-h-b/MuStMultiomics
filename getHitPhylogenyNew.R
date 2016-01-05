# this script is used in the MuSt project to find the taxonomic annotation and mOTU of genes (of one COG) found in a BLAST search

#inputs: as argument in the call - COG, the COG ID of the genes of interest
#       uses - the mapping files "mOTU.v1.padded.motu.linkage.map", "mOTU.v1.map.txt" and "mOTU-LG.v1.annotations.txt" of the mOTU database and
#            - a prepared file with the best hit output of BLASTN "<COG>.bestHits.tsv"  
#output: a table with the query gene and the taxonomy of the best hit (written to a file named "<COG>.bestHitPhylogeny.tsv")
# 

#written by Anna Heintz-Buschart (March 2015)

args<-commandArgs(TRUE)
COG <- args[1]

phyl <- read.delim("../../mOTU.v1.padded.motu.linkage.map",stringsAsFactors=F)
map <- read.delim("../../mOTU.v1.map.txt",stringsAsFactors=F)
mapPhyl <- read.delim("../../mOTU-LG.v1.annotations.txt",stringsAsFactors=F)

bh <- read.delim(paste(COG,".bestHits.tsv",sep=""),header=F,comment.char="#",stringsAsFactors=F)
colnames(bh) <- c("queryID", "subjectID", "percIdentity", "alignmentLength", "mismatches", "gapOpens", "q.start", "q.end", "s.start", "s.end", "evalue", "bitScore")

bh$subjectID <- sapply(bh$subjectID,function(x) paste(unlist(strsplit(x,split=".",fixed=T))[-c(length(unlist(strsplit(x,split=".",fixed=T)))-1:0)],sep=".",collapse="."))
bh <- merge(bh,map,by.x=2,by.y=1,all.x=T)
bh <- merge(bh,phyl[,c(1,11)],by.x="mOTU",by.y=1,all.x=T)
bh <- merge(bh,mapPhyl,by.x=ncol(bh),by.y=1,all.x=T)
write.table(bh,paste(COG,".bestHitPhylogeny.tsv",sep=""),quote=F,sep="\t",row.names=F)
