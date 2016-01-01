# this script parses organism_RefSeq.tab files from MG-RAST, to get taxonomic annotations for genes
# all annotations with scores within 80% of the top score are used (no lower threshold is applied at this step)
# the script relies on the NCBI taxdump and parse_taxbrowser_MGRAST.py, which should be present in the taxdump directory
# the MG-RAST data needs to be present in one directory, and a table with the MG-RAST IDs and the project IDs needs to be prepared beforehand
# SET ALL THE PATHS in lines 9-11

# Anna Heintz-Buschart, version of November 2015

wd <- "refseqdirectory/" # # # set your working directory containing all the organism_RefSeq.tab files from MG-RAST
taxdir <- "taxdump/" # # # set the pathway to the NCBI taxdump here (with a slash in the end)
idsFile <- "MGRastIDs.txt" # # # table with a translation between the IDs from MG-RAST (used in the file name of the .tab file) and your prefered IDs

setwd(wd)

extrContig <- function(contigNames,RASTid) gsub("_[0-9]+_[0-9]+_.$","",gsub(paste("mgm",RASTid,".3.",sep=""),"",contigNames))

orgAnnoTPF <- function(RASTid,sample,taxdir){
  annoTabO <- read.delim(paste("mgm",RASTid,".3_organism_RefSeq.tab",sep=""),stringsAsFactors=F)
  if(length(grep("Download complete",annoTabO[nrow(annoTabO),]))==0) stop(paste("The file", RASTid, ".3_organism_RefSeq.tab is incomplete."))
  annoTabO <- annoTabO[-nrow(annoTabO),]
  if(nrow(annoTabO)==0) stop(paste("The file", RASTid, ".3_organism_RefSeq.tab is empty."))
    
  annoTabO[,1] <- extrContig(annoTabO[,1],RASTid)
  annoTabOc <- annoTabO[,-2] 
  rm(annoTabO)
  annoTabOc[ is.na(annoTabOc$bit.score)] <- 0
  
  agCon <- aggregate(annoTabOc[,c(11,12)],list(annoTabOc$query.sequence.id),function(x)x)
  strain <- apply(agCon[,2:3],1,function(x) paste(unique(unlist(strsplit((x[2][[1]][x[1][[1]]>=0.8*max(x[1][[1]])]),split=";"))),
                                                  sep=";",collapse=";") )
  
  contigs <- unique(annoTabOc$query.sequence.id)
  contigs <- contigs[order(contigs)]
  
  oTabS <- data.frame("contig"=contigs,"strain"=strain,stringsAsFactors=F)
  oTabS[,ncol(oTabS)+1] <- sapply(oTabS$strain,function(x) length(unlist(strsplit(x,split=";"))))
  colnames(oTabS)[ncol(oTabS)] <- "hits.number"
  
  singOrg <- oTabS[ oTabS$hits.number ==1, ][,1:2]
  multOrg <- oTabS[ oTabS$hits.number !=1, ][,1:2]
  
  write.table(multOrg[,2],paste(RASTid,"multOut.txt",sep=""),row.names=F,col.names=F,quote=F)
  write.table(singOrg[,2],paste(RASTid,"singOut.txt",sep=""),row.names=F,col.names=F,quote=F)
  system(paste("python ",taxdir,"parse_taxbrowser_MGRAST.py ",taxdir, " ",wd,"/", RASTid,"multOut.txt ",wd,"/",RASTid,"singOut.txt ",wd,"/",RASTid,"multsingPhylogeny.txt",sep=""))
  
  singedOrg <- rbind(multOrg,singOrg)
  singedOrgPhyl <- data.frame("gene"=singedOrg[,1],read.delim(paste(RASTid,"multsingPhylogeny.txt",sep=""),stringsAsFactors=F),stringsAsFactors=F)
  
  write.table(singedOrgPhyl[singedOrgPhyl$kingdom=="Eukaryota"&singedOrgPhyl$order!="Primates",],paste("eukaryota/",sample,".eukaryota.tsv",sep=""),sep="\t",quote=F,row.names=F)
  write.table(singedOrgPhyl[singedOrgPhyl$kingdom=="Viruses",],paste("virus/",sample,".viruses.tsv",sep=""),sep="\t",quote=F,row.names=F)
  return()
}


Mall <- read.delim(idsFile,stringsAsFactors=F,header=F) 

Mall[,3] <- gsub(".3","",Mall[,3],fixed=T)

for(i in 1:nrow(Mall)){
  orgAnnoTPF(Mall[i,3],gsub("M0","M",Mall[i,2]),taxdir)
}

