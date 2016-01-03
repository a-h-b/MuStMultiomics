#this script gets all genes with a function of interest as best annotation 
#### and creates a number of plots which show from which genomes this function is expressed in the different samples and on the different omic levels
#### it also returns a workspace with the data

# it takes 3 ARGUMENTs when called: the function of interest, which mOTU annotation (best hit out of the mOTUs found at reads level ("mOTUpresent") or
#### of all ("mOTUbest")) to use and whether to order the plots by whether the donors of the samples have T1DM ("T1DM") or belong to a group defined in another file ("BG")
# the name of the function of interest is used to create a directory which houses all the output plots

# the script is constructed in two parts: the part that accesses the database and the part that makes the plots. The plotting needs a lot of additional informations
### and also the script 140510_heatmap2.R. The database access is found in lines 14-15, and 34-76. The rest is plotting.

# written by Anna Heintz-Buschart, this version is from October 2015

args<-commandArgs(TRUE)
funOI <- args[1]
####
motuCh <- args[2]
ord <- args[3]

source("140510_heatmap2.R")
clusterInfo <- readRDS("allClusterInfo.RDS")
motu <- readRDS("motuAnnotation.RDS")
bigGroup <- readRDS("bigGroup.RDS")
mapStats <- readRDS("mappedReads.RDS")
protStats <- read.delim("allStats.tsv",stringsAsFactors=F,row.names=1)
colnames(protStats) <- gsub(".V","-V",colnames(protStats),fixed=T)
miMetaA <- readRDS("miMetaCombi.RDS")

library(gplots)
library(RColorBrewer)
library(vegan)
########

library(rmongodb)

getExprTab <- function(funOI){
  if(mongo.is.connected(mongo)) {
    unwindGenes <- mongo.bson.from.JSON('{"$unwind": "$genes"}')
    matchAnno <- mongo.bson.from.JSON(paste('{"$match": {"genes.bestAnnotation": "',funOI,'"} }',sep=""))
    matchProt <- mongo.bson.from.JSON(paste('{"$match": {"genes.proteinIdentification": "uniquely"} }',sep=""))
    groupDNARNAname <- mongo.bson.from.JSON('{"$group": {"_id": "$_id", "aveCovDNA": {"$push": "$genes.aveCovDNA"},"sample":{"$push":"$sample"},"gene":{"$push":"$genes.gene"},"cluster":{"$push": "$cluster"},"aveCovRNA":{"$push": "$genes.aveCovRNAfw"}}}')
    projectDNARNAname <- mongo.bson.from.JSON('{"$project": {"_id": 0, "gene":1,"aveCovDNA": 1, "cluster" :1,"aveCovRNA":1,"sample":1}}')
    groupProtname <- mongo.bson.from.JSON('{"$group": {"_id": "$_id", "gene": {"$push": "$genes.gene"},"sample":{"$push":"$sample"},
                                        "proteinIdentification":{"$push":"$genes.proteinIdentification"},
                                        "proteinArea":{"$push":"$genes.proteinArea"}}}') 
    projectProtname <- mongo.bson.from.JSON('{"$project": {"_id": 0, "gene":1,"proteinIdentification":1,"proteinArea":1,"sample":1}}')
    genesA <- mongo.bson.to.list(mongo.aggregation(mongo,coll,list(matchAnno,unwindGenes,matchAnno,groupDNARNAname,projectDNARNAname)))$result
    genesP <- mongo.bson.to.list(mongo.aggregation(mongo,coll,list(matchAnno,unwindGenes,matchAnno,matchProt,groupProtname,projectProtname)))$result
    if(length(genesA)==0){ 
      warning("No genes found")
      return(NULL)
    }
    res <- do.call(rbind,lapply(genesA,function(x) cbind(x$sample,x$gene,x$cluster,x$aveCovDNA,x$aveCovRNA)))
    res <- data.frame("sample"=res[,1],"gene"=res[,2],"cluster"=res[,3],"aveCovDNA"=as.numeric(res[,4]),"aveCovRNA"=as.numeric(res[,5]),
                      stringsAsFactors=F)
    if(length(genesP)>0){
      pFeat <- do.call(rbind,lapply(genesP,function(x) cbind(x$sample,x$gene, x$proteinArea)))
      pFeat[sapply(pFeat[,3],length)>1,3] <- sapply(pFeat[sapply(pFeat[,3],length)>1,3], mean)
      pFeat <- data.frame("sample"=unlist(pFeat[,1]),"gene"=unlist(pFeat[,2]),"proteinArea"=as.numeric(unlist(pFeat[,3])),stringsAsFactors=F)
      aFeat <- merge(res,pFeat,by=c(1,2),all.x=T)
    }else{
      aFeat <- res
      aFeat$proteinArea <- 0
    }
    aFeat[is.na(aFeat)] <- 0
    return(aFeat)
  } else {
    stop("Mongo is not connected.")
  }
}

mongo <- mongo.create()
db <- "mydb" 
coll <- "mydb.must" 
exTab <- getExprTab(funOI)
mongo.destroy(mongo)

#######

mapFac <- mapStats[,2:3]
rownames(mapFac) <- mapStats[,1]
mapFac$DNAreadsOnContigs <- 1/(mapStats[,2]/mean(mapStats[,2]))
mapFac$RNAreadsFwOnGenes <- 1/(mapStats[,3]/mean(mapStats[,3]))
miMeta <- miMetaA[-10,]
McolA <- c(brewer.pal(7,"Blues")[3:7],brewer.pal(7,"Oranges")[3:7],brewer.pal(6,"Purples")[3:6],brewer.pal(8,"Greens")[3:8]) #all 4 families
Mcol <- McolA[-c(5,11)]
visFacA <- c(rep(1:4,each=3),rep(5,times=2),rep(6:9,each=3),rep(10:11,each=2),rep(12:13,each=3),rep(14:15,each=2),rep(16:20,each=3)) #all 4 families
visFac <- c(rep(1:3,each=3),rep(4,times=2),rep(5,times=3),rep(6,times=2),rep(7:8,each=3),rep(9,times=2),10:13,rep(14,times=2),15,rep(16,times=2),17,rep(18,times=2))
rownames(miMeta) <- gsub("-0","-",gsub("M-","M",rownames(miMeta)))
indiv <- c(paste("M01.",c(1:4),sep=""),paste("M02.",c(1:5),sep=""),paste("M03.",c(3:5),sep=""),
           paste("M04.",c(1:6),sep=""))
hmgrey <- colorRampPalette(brewer.pal(9,"Greys"),bias=2.5)(256)
DNAhm <- colorRampPalette(c("white",rgb(0,204,204,maxColorValue=255),rgb(0,102,102,maxColorValue=255)))(256)
RNAhm <- colorRampPalette(c("white",rgb(204,0,204,maxColorValue=255),rgb(102,0,102,maxColorValue=255)))(256)
Prothm <- colorRampPalette(c("white",rgb(204,204,0,maxColorValue=255),rgb(102,102,0,maxColorValue=255)))(256)
rathm <- colorRampPalette(c(rgb(0,102,102,maxColorValue=255),rgb(0,204,204,maxColorValue=255),"white",rgb(204,0,204,maxColorValue=255),rgb(102,0,102,maxColorValue=255)))(256)
plotExprClus <- function(exprTab,cI=clusterInfo,oOI=funOI,motuChoice=motuCh,retVal="plot"){
  require(RColorBrewer)
  for(lib in sort(unique(exprTab$sample))){
    popRNACov <- aggregate(exprTab$aveCovRNA[exprTab$sample==lib],list(exprTab$cluster[exprTab$sample==lib]),sum)
    colnames(popRNACov) <- c("cluster","cumCovRNA")
    popGenes <- aggregate(exprTab$aveCovRNA[exprTab$sample==lib],list(exprTab$cluster[exprTab$sample==lib]),length)
    colnames(popGenes) <- c("cluster","geneNo")
    pops <- merge(popRNACov,cI[cI$sample==lib,c("cluster","aveCov","uniqueEss",motuChoice)],by=1,all.x=T)
    pops$motuUnanimous <- ifelse(grepl(";",pops[[motuChoice]])|pops[[motuChoice]]=="","uncertain",
                                 sapply(pops[[motuChoice]],function(x) unlist(strsplit(x,split="\\("))[1]))
    #pops$motuUnanimous[pops$motuUnanimous==""] <- ""
    if(any(!pops$cluster %in% c("N","S"))){
      exTab <- pops[!pops$cluster %in% c("N","S"),]
      exTab <- merge(exTab,popGenes,by="cluster")
      exTab <- exTab[order(exTab$cumCovRNA,decreasing=T),]
      exTab$cluster <- paste(exTab$cluster," (",exTab$geneNo,ifelse(exTab$geneNo==1," gene)"," genes)"),"\n",ifelse(exTab$motuUnanimous %in% c("uncertain",""),exTab$motuUnanimous,gsub("SpeciesCluster of ","",sapply(exTab$motuUnanimous,function(x)motu$SpeciesCluster[motu$ID==x]))),sep="")
      exTab <- exTab[!is.na(exTab$cumCovRNA),]
      if(any(pops$cluster %in% c("N","S"))){
        exTab <- rbind(exTab[,1:4],c(paste(sum(popGenes$geneNo[popGenes$cluster %in% c("N","S")]),
                                           ifelse(sum(popGenes$geneNo[popGenes$cluster %in% c("N","S")])==1,"other gene","other genes")),
                                     sum(popRNACov$cumCovRNA[popRNACov$cluster %in% c("N","S")]),1))
      }
    }else{
      exTab <- data.frame("cluster"=paste(sum(popGenes$geneNo[popGenes$cluster %in% c("N","S")]),
                                          ifelse(sum(popGenes$geneNo[popGenes$cluster %in% c("N","S")])==1,"other gene","other genes")),
                          "cumCovRNA"=sum(popRNACov$cumCovRNA[popRNACov$cluster %in% c("N","S")]),"aveCov"=NA,"uniqueEss"=0)
    }
    if("plot" %in% retVal){
      maxy <- 1.1*max(as.numeric(c(exTab$cumCovRNA)))
      par(mar=c(3,10,0.5,0.5),mgp=c(1.9,0.6,0))
      barplot(as.numeric(exTab$cumCovRNA),names.arg=exTab$cluster,las=2,xlim=c(0,maxy),horiz=T,cex.names=0.6,
              xlab=paste("metaT coverage depth of",oOI),
              col=colorRampPalette(brewer.pal(11,"Spectral"))(109)[109:1][as.numeric(exTab$uniqueEss)+1],
              cex.axis=0.8)
      mtext(paste("populations in sample",lib),2,8.7,cex=1.2)
      if(any(popGenes$cluster %in% c("N","S"))){
        popCount <- nrow(exTab)-1
        geneCount <- popGenes$geneNo[popGenes$cluster %in% c("N","S")]
        barplot(cbind(matrix(0,nrow=sum(geneCount),ncol=popCount),
                      sort(exprTab$aveCovRNA[exprTab$sample==lib&exprTab$cluster %in% c("N","S")],decreasing=T)),
                names.arg=rep("",popCount+1),axes=F,add=T,horiz=T)
      }
      if(any(!grepl("other",exTab$cluster))){
        maxx <- 1.1*max(as.numeric(exTab$aveCov[!grepl("other",exTab$cluster)]))
        maxy <- 1.1*max(as.numeric(exTab$cumCovRNA[!grepl("other",exTab$cluster)]))
        par(mar=c(3,3,1.1,0.5))
        plot(as.numeric(exTab$aveCov[!grepl("other",exTab$cluster)]),as.numeric(exTab$cumCovRNA[!grepl("other",exTab$cluster)]),
             las=1,xlim=c(0,maxx),ylim=c(0,maxy),xlab="metaG coverage depth of cluster",ylab=paste("metaT coverage depth of",oOI),
             col=colorRampPalette(brewer.pal(11,"Spectral"))(109)[109:1][as.numeric(exTab$uniqueEss[!grepl("other",exTab$cluster)])+1],
             cex.axis=0.8,cex=sqrt(as.numeric(exTab$geneNo[!grepl("other",exTab$cluster)])),pch=16)
        mtext(lib,3,0,cex=1.2)
        text(as.numeric(exTab$aveCov[!grepl("other",exTab$cluster)]),as.numeric(exTab$cumCovRNA[!grepl("other",exTab$cluster)]),
             labels=exTab$cluster[!grepl("other",exTab$cluster)],adj=c(-0.2,-0.2),cex=0.8,
             col=colorRampPalette(brewer.pal(11,"Spectral"))(109)[109:1][as.numeric(exTab$uniqueEss[!grepl("other",exTab$cluster)])+1])
      }
    }
  }
  if("max" %in% retVal) retList <- exprTab$cluster[which.max(exprTab$aveCovRNA)]
  if("exprTabMotu" %in% retVal) retList <- merge(exprTab,cI[,c("sample","cluster",motuChoice)],by=c("sample","cluster"))
  return(retList)
}


outDir <- paste("./",funOI,sep="")
dir.create(outDir)
pdf(paste(outDir,"/metaTcovAllClusters.pdf",sep=""),width=3.5,height=3.5,pointsize=8)
resL <- plotExprClus(exTab,retVal=c("exprTabMotu","plot"))
dev.off()
if(!is.null(resL)){
  resL$motuUnanimous <- ifelse(grepl(";",resL$motuPresent)|resL$motuPresent == "","uncertain",
                               sapply(resL$motuPresent,function(x) unlist(strsplit(x,split="\\("))[1]))
  resL$motuUnanimous[is.na(resL$motuUnanimous)] <- "uncertain"
  
  resLclus <- tapply(resL$cluster,list(resL$motuUnanimous,resL$sample),function(x) length(unique(x)))
  resLclus[is.na(resLclus)] <- 0
  resLgene <- tapply(resL$aveCovDNA,list(resL$motuUnanimous,resL$sample),length)
  resLgene[is.na(resLgene)] <- 0
  
  mapFac <- mapFac[rownames(mapFac) %in% unique(resL$sample),]
  visFac <- visFac[sort(unique(clusterInfo$sample)) %in% unique(resL$sample)]
  
  resLDNA <- tapply(resL$aveCovDNA,list(resL$motuUnanimous,resL$sample),sum)
  resLDNA[is.na(resLDNA)] <- 0
  resLDNA <- t(apply(resLDNA,1,function(x) x*mapFac$DNAreadsOnContigs))
  
  resLRNA <- tapply(resL$aveCovRNA,list(resL$motuUnanimous,resL$sample),sum)
  resLRNA[is.na(resLRNA)] <- 0
  resLRNA <- t(apply(resLRNA,1,function(x) x*mapFac$RNAreadsFwOnGenes))
  
  resLProtein <- tapply(resL$protein,list(resL$motuUnanimous,resL$sample),sum)
  resLProtein[is.na(resLProtein)] <- 0
  protFac <- as.vector(t(protStats[rownames(protStats)=="totalArea",colnames(protStats) %in% resL$sample]))
  resLProtein <- t(apply(resLProtein,1,function(x) x/protFac))
  if(nrow(resLgene)>1){
    resLgene <- resLgene[order(rowSums(resLRNA),decreasing=T),]
    resLclus <- resLclus[order(rowSums(resLRNA),decreasing=T),]
    resLDNA <- resLDNA[order(rowSums(resLRNA),decreasing=T),]
    resLProtein <- resLProtein[order(rowSums(resLRNA),decreasing=T),]
    resLRNA <- resLRNA[order(rowSums(resLRNA),decreasing=T),]
    

    trsl <- sapply(rownames(resLgene),function(x) if(x=="uncertain") x else(motu$SpeciesCluster[motu$ID == x]))
    trsl <- gsub("SpeciesCluster of ","",trsl)
    
    pdf(paste(outDir,"/heatmaps_all.pdf",sep=""),height=4.9,width=7.8,pointsize=8)
    if(any(resLgene>0)){
      if(ord=="T1DM"){
        ordervec <- c(which(miMeta$DIABETESTY1[visFac]=="Yes"),which(miMeta$DIABETESTY1[visFac]!="Yes"))
      }else if(ord=="BG"){
        ordervec <- c(which(bigGroup[visFac]==1),which(bigGroup[visFac]==2))
      }else{
        ordervec <- c(1:ncol(resLgene))
      }
      heatmap.2a(resLgene[,ordervec],keyName="number of genes",
                 trace="n",Colv="none",col=hmgrey,Rowv="none",density.info="n",labRow=trsl,dendrogram="none",margins=c(4,25),
                 ColSideColors=rbind(Mcol[visFac],c("black","white")[1+as.numeric(miMeta$DIABETESTY1[visFac]=="No")])[,ordervec])
    }
    if(any(resLclus>0)){
      heatmap.2a(resLclus[,ordervec],keyName="number of clusters",
                 trace="n",Colv="none",col=hmgrey,Rowv="none",density.info="n",labRow=trsl,dendrogram="none",margins=c(4,25),
                 ColSideColors=rbind(Mcol[visFac],c("black","white")[1+as.numeric(miMeta$DIABETESTY1[visFac]=="No")])[,ordervec])
    }
    if(any(resLDNA>0)){
      heatmap.2a(resLDNA[,ordervec],keyName="cluster cov metaG",
                 trace="n",Colv="none",col=DNAhm,Rowv="none",density.info="n",labRow=trsl,dendrogram="none",margins=c(4,25),
                 ColSideColors=rbind(Mcol[visFac],c("black","white")[1+as.numeric(miMeta$DIABETESTY1[visFac]=="No")])[,ordervec])
    }
    if(any(resLRNA>0)){
      heatmap.2a(resLRNA[,ordervec],keyName="gene cov metaT",
                 trace="n",Colv="none",col=RNAhm,Rowv="none",density.info="n",labRow=trsl,dendrogram="none",margins=c(4,25),
                 ColSideColors=rbind(Mcol[visFac],c("black","white")[1+as.numeric(miMeta$DIABETESTY1[visFac]=="No")])[,ordervec])
    }
    if(any(resLRNA>0)&any(resLDNA>0)){
      resRat <- log10(resLRNA[,ordervec]+0.001)-log10(resLDNA[,ordervec]+0.001)
      heatmap.2a(resRat[,ordervec],keyName="log10 metaT/metaG",symbreaks=T,
                 trace="n",Colv="none",col=rathm,Rowv="none",density.info="n",labRow=trsl,dendrogram="none",margins=c(4,25),
                 ColSideColors=rbind(Mcol[visFac],c("black","white")[1+as.numeric(miMeta$DIABETESTY1[visFac]=="No")])[,ordervec])
    }
    if(any(resLProtein>0)){
      heatmap.2a(resLProtein[,ordervec],keyName="protein rel quant",
                 trace="n",Colv="none",col=Prothm,Rowv="none",density.info="n",labRow=trsl,dendrogram="none",margins=c(4,25),
                 ColSideColors=rbind(Mcol[visFac],c("black","white")[1+as.numeric(miMeta$DIABETESTY1[visFac]=="No")])[,ordervec])
    }
    dev.off()
    
    top <- apply(resLRNA,1,function(x) max(x/colSums(resLRNA),na.rm=T))>0.1
    
    pdf(paste(outDir,"/heatmaps_10perc.pdf",sep=""),height=4.9,width=7.8,pointsize=8)
    if(any(resLgene[top,]>0)){
      heatmap.2a(resLgene[top,ordervec],keyName="number of genes",
                 trace="n",Colv="none",col=hmgrey,Rowv="none",density.info="n",labRow=trsl[top],dendrogram="none",margins=c(4,25),
                 ColSideColors=rbind(Mcol[visFac],c("black","white")[1+as.numeric(miMeta$DIABETESTY1[visFac]=="No")])[,ordervec])
    }
    if(any(resLclus[top,]>0)){
      heatmap.2a(resLclus[top,ordervec],keyName="number of clusters",
                 trace="n",Colv="none",col=hmgrey,Rowv="none",density.info="n",labRow=trsl[top],dendrogram="none",margins=c(4,25),
                 ColSideColors=rbind(Mcol[visFac],c("black","white")[1+as.numeric(miMeta$DIABETESTY1[visFac]=="No")])[,ordervec])
    }
    if(any(resLDNA[top,]>0)){
      heatmap.2a(resLDNA[top,ordervec],keyName="cluster cov metaG",
                 trace="n",Colv="none",col=DNAhm,Rowv="none",density.info="n",labRow=trsl[top],dendrogram="none",margins=c(4,25),
                 ColSideColors=rbind(Mcol[visFac],c("black","white")[1+as.numeric(miMeta$DIABETESTY1[visFac]=="No")])[,ordervec])
    }
    if(any(resLRNA[top,]>0)){
      heatmap.2a(resLRNA[top,ordervec],keyName="gene cov metaT",
                 trace="n",Colv="none",col=RNAhm,Rowv="none",density.info="n",labRow=trsl[top],dendrogram="none",margins=c(4,25),
                 ColSideColors=rbind(Mcol[visFac],c("black","white")[1+as.numeric(miMeta$DIABETESTY1[visFac]=="No")])[,ordervec])
    }
    if(any(resLRNA[top,]>0)&any(resLDNA[top,]>0)){
      resRat <- log10(resLRNA[top,ordervec]+0.001)-log10(resLDNA[top,ordervec]+0.001)
      heatmap.2a(resRat[,ordervec],keyName="log10 metaT/metaG",symbreaks=T,
                 trace="n",Colv="none",col=rathm,Rowv="none",density.info="n",labRow=trsl,dendrogram="none",margins=c(4,25),
                 ColSideColors=rbind(Mcol[visFac],c("black","white")[1+as.numeric(miMeta$DIABETESTY1[visFac]=="No")])[,ordervec])
    }
    if(any(resLProtein[top,]>0)){
      heatmap.2a(resLProtein[top,ordervec],keyName="protein rel quant",
                 trace="n",Colv="none",col=Prothm,Rowv="none",density.info="n",labRow=trsl[top],dendrogram="none",margins=c(4,25),
                 ColSideColors=rbind(Mcol[visFac],c("black","white")[1+as.numeric(miMeta$DIABETESTY1[visFac]=="No")])[,ordervec])
    }
    dev.off()
    save.image(paste(outDir,"/WS.Rdata",sep=""))
  }
  
}
