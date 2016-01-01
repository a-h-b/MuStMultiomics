# R script to reconstruct metabolic network based on KOs detected in a dataset
# needs: the library BioNet; the analysis function in this script uses FastHeinz within R

# input: several tables obtained from rest.kegg.jp (see comments below)
# outputs: a network as graphNEL R-object (koGraph.RDS) and as comma separated edge table (NW.csv; can be imported into Cytoscape for example)
#           and a file with all KOs in the network (KOs_in_NW.tsv) and all the edges (allEdges.tsv)
# the analysis function can also create a pdf file with a plot of the highest scoring subnetwork from BioNet analysis

# written by Anna Heintz-Buschart, based on an idea by Hugo Roume
# adapted from code used in the publication "Comparative integrated omics: identification of key functionalities in microbial community-wide metabolic networks" by Roume, H., Heintz-Buschart, et al. (2015). npj Biofilms and Microbiomes, 1, 15007. doi:10.1038/npjbiofilms.2015.7

# this is version 1.1 from June 2015 (commenting changed in Dec 2015)

library(BioNet)
substrRight <- function(x, n){
  substr(x, nchar(x)-n+1, nchar(x))
}

#names of KOs - to be used later 
kon <- read.delim("ko2des_clean.txt",stringsAsFactors=F,header=F,quote="")
colnames(kon) <- c("koID","description")


#which pathway does a ko belong to?, keep only metabolic pathways
pw <- read.delim("150630_kegg_ko2pathway.txt",header=F,stringsAsFactors=F) #this list comes from http://rest.kegg.jp/link/Ko/pathway
pw <- pw[c(grep("path:ko00",pw[,1]),grep("path:ko01",pw[,1])),]
pw[,1] <- gsub("path:","",pw[,1])
pw[,2] <- gsub("ko:","",pw[,2])
colnames(pw) <- c("pw","ko")

#read ko mapping to reaction
rn2ko <- read.delim("150630_kegg_RN2KO.txt",header=F,stringsAsFactors=F) #this list comes from http://rest.kegg.jp/link/Ko/rn
rn2ko[,1] <- gsub("rn:","",rn2ko[,1])
rn2ko[,2] <- gsub("ko:","",rn2ko[,2])
colnames(rn2ko) <- c("rn","ko")

koOI <- unique(pw$ko)

#merge kos from metabolic pathways with reactions
rnkoOI <- merge(koOI,rn2ko,by.x=1,by.y="ko")
rnkoOI[,1] <- as.character(rnkoOI[,1])
colnames(rnkoOI)[1] <- "ko"

#read mapping of reaction pairs to reactions
r2rp <- read.delim("150630_kegg_RP2RN.txt",header=F,stringsAsFactors=F) #this list comes from http://rest.kegg.jp/link/rn/rp
r2rp[,1] <- gsub("rp:","",r2rp[,1])
r2rp[,2] <- gsub("rn:","",r2rp[,2])
colnames(r2rp) <- c("rp","rn")

#merge reaction pairs with kos via reactions
rprnkoOI <- merge(rnkoOI,r2rp,by.x=2,by.y=2)

#read mapping of reaction pairs to reaction class
rp2rc <- read.delim("150630_kegg_rc2rp.txt",header=F,stringsAsFactors=F) #this list comes from http://rest.kegg.jp/link/rp/rc/
rp2rc[,1] <- gsub("rc:","",rp2rc[,1])
rp2rc[,2] <- gsub("rp:","",rp2rc[,2])
colnames(rp2rc) <- c("rc","rp")

#merge reaction class with kos via reaction pairs
rpcrnkoOI <- merge(rprnkoOI,rp2rc,by.x=3,by.y=2)

#merge reaction class with kos without reaction pairs
rn2rc <- read.delim("150630_kegg_rc2rn.txt",header=F,stringsAsFactors=F) #this list comes from http://rest.kegg.jp/link/rc/rn
rn2rc[,1] <- gsub("rn:","",rn2rc[,1])
rn2rc[,2] <- gsub("rc:","",rn2rc[,2])
colnames(rn2rc) <- c("rn","rc")

crnkoOI <- merge(rnkoOI,rn2rc,by.x=2,by.y=1)

kornrpOI <- merge(crnkoOI,rpcrnkoOI,by.x=c(1,2,3),by.y=c(2,3,4))
kornrpOI <- unique(kornrpOI)

#read mapping of reaction pairs to compounds
rp2c <- read.delim("150630_kegg_RP2cpd.txt",header=F,stringsAsFactors=F) #this list comes from http://rest.kegg.jp/list/rp
rp2cT <- rp2c
rp2cT[,1] <- gsub("rp:","",rp2cT[,1])
#split compounds
rp2cT[,2] <- sapply(rp2cT[,2],function(x){a <- strsplit(x,split="_"); unlist(a)})[1,]
rp2cT[,3] <- sapply(rp2c[,2],function(x){a <- strsplit(x,split="_"); unlist(a)})[2,]
colnames(rp2cT) <- c("rp","cpd1","cpd2")

#merge kos with compounds via reactions and reaction pairs - keep only ko and compounds
crprnkoOI <- merge(kornrpOI,rp2cT,by.x=4,by.y=1)
ckoOI <- crprnkoOI[,c(3,5,6)]
ckoOI <- unique(ckoOI,MARGIN=1)
ckoOIc2 <- data.frame("ko"= rep(ckoOI$ko,2),"C"=c(ckoOI$cpd1,ckoOI$cpd2),stringsAsFactors=F)
ckoOIc2 <- unique(ckoOIc2,MARGIN=1)
# no more water in the network? ckoOIc2[ckoOIc2$C=="C00001",]


#list all compounds of a ko
lckoOI <- aggregate(ckoOIc2$C,list(ckoOIc2$ko),paste)
lck <- as.vector(lckoOI[,2],mode="list")
names(lck) <- lckoOI[,1]
#-> 3953 KOs

#should not be necessary:
for(i in 1:length(lck)){
  lck[[i]] <- unique(lck[[i]])
}

lcko <- lck
#combine KOs with the same compounds, but only if the KO numbers are adjacent, because then they are subunits of the same enzyme!
i <- 1
while(i <= length(lck)){
  j <- i+1
  while(j <= length(lck)){
    if(setequal(lck[[i]],lck[[j]]) & (as.numeric(substrRight(names(lck)[i],5))+1==as.numeric(substrRight(names(lck)[j],5)))){
      names(lck)[i] <- paste(names(lck)[i],names(lck)[j],sep= "-")
      lck <- lck[-j]
    }else j <- j+1
  }
  i <- i+1
}
#->  3367 combined KOs


### join KOs with edges (include which compound)
n1L <- vector()
n2L <- vector()
edgesL <- list()
allEdL <- vector()
eCoL <- vector()
gMat <- matrix(0,nrow=length(lck),ncol=length(lck),dimnames=list(names(lck),names(lck)))
for(i in 1:(length(lck)-1)){
  for(j in (i+1):length(lck)){
    ise <- intersect(lck[[i]],lck[[j]])
    noCo <- length(ise)
    if(noCo>0){
      gMat[i,j] <- 1
      gMat[j,i] <- 1
      n1L <- append(n1L,names(lck)[i])
      n2L <- append(n2L,names(lck)[j])
      edgesL[[length(edgesL)+1]] <- ise
      allEdL <- append(allEdL,ise)
      eCoL <- append(eCoL,noCo)
    }
  }
}
# length(unique(c(n1L,n2L))) gives number of nodes
#-> 3353 connected, 14 unconnected

### calculate number of occurences for every compound that makes an edge
allEdNoL <- table(allEdL)

#read names of compounds
cname <- read.delim("150630_kegg_cpdNames.txt",header=F,stringsAsFactors=F) #this list comes from http://rest.kegg.jp/list/cpd
cname[,1] <- gsub("cpd:","",cname[,1])
allEdTab <- data.frame("cID"=names(allEdNoL)[order(allEdNoL,decreasing=T)],"abundance"=allEdNoL[order(allEdNoL,decreasing=T)])
allEdTab <- merge(allEdTab,cname,by.x=1,by.y=1,sort=F)
write.table(allEdTab,"allEdges.tsv",sep="\t",row.names=F)

### edge names
edgeNamesL <- vector()
for(i in 1:length(n1L)){
  if(length(edgesL[[i]])==1){
    edgeNamesL <- append(edgeNamesL,edgesL[[i]])
  }else{
    helpEdge <-  edgesL[[i]][1]  
    for(k in 2:length(edgesL[[i]])){
      helpEdge <- paste(helpEdge,edgesL[[i]][k],sep="-")
    }
    edgeNamesL <- append(edgeNamesL,helpEdge)
  }
}

### combine into data frame
allEdgesNodesL <- data.frame("node1"=n1L,"node2"=n2L,"edgename"=edgeNamesL,"compoundNo"=eCoL)
write.table(allEdgesNodesL,"NW.csv",sep=",",row.names=F,quote=F)

### Use BioNet ####
koGraphAM <- graphAM(adjMat=gMat,edgemode="undirected")
koGraph <- as(koGraphAM,"graphNEL")
saveRDS(koGraph,"koGraph.RDS")
write.table(nodes(koGraph),"KOs_in_NW.tsv",row.names=F,col.names=F,sep="\t")


### to get a network with just nodes represented by detected KOs (eg the rownames of kAbundNT):
# kAbundNT <- readRDS("kAbundNT.RDS")
# potVert <- nodes(koGraph)[which(sapply(nodes(koGraph),function(x)length(intersect(unique(unlist(strsplit(x,split="-"))),rownames(kAbundNT)))>0))]
# 
# mustGraphMain <- largestComp(subNetwork(potVert,koGraph))
# mustIGraph <- igraph.from.graphNEL(mustGraphMain,name=T,weight=T)
# saveRDS(mustGraphMain,"mustGraph.RDS")
# saveRDS(mustIGraph,"mustIGraph.RDS")


### function to load data and plot module and return nodes;
## it can be used on tables containing node names, p-values (column name pvalue) and logFC (column name log2FoldChange), like the output of DESeq2
annaMustPlotNW <- function(file,outname="module.pdf",fdr=0.05,graph=koGraph,w2d=F){
  data <- read.delim(file,stringsAsFactors=F,row.names=1)
  subnet <- subNetwork(rownames(data)[!is.na(data$pvalue)],graph)
  subnet <- rmSelfLoops(subnet)
  pval <- data$pvalue[!is.na(data$pvalue)]
  names(pval) <- rownames(data)[!is.na(data$pvalue)]
  fb <- fitBumModel(pval, plot = FALSE)
  scores <- scoreNodes(subnet, fb, fdr = fdr)
  module <- runFastHeinz(subnet, scores)
  logFC <- data$log2FoldChange[!is.na(data$pvalue)]
  names(logFC) <- rownames(data)[!is.na(data$pvalue)]
  if(w2d){
    pdf(outname,width=3.5,height=3.5,pointsize=8)
    plotModule(module, scores = scores, diff.expr = logFC)
    dev.off()
  }
  plotModule(module, scores = scores, diff.expr = logFC)
  return(nodes(module))
}
