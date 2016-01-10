# this script can be used to retrieve metagenomic, metatranscriptomic and metaproteomic contributions of each reconstructed population-level genomes
### (and all contigs >= 1000 bp) to all nodes in the reconstructed network represented by "150705_KOs_in_NW.tsv" from the Mongo DB
#
#inputs: in addition to "150705_KOs_in_NW.tsv", this script needs an R-object with information on all reconstructed population-level genomes, "allClusterInfo.RDS"
# outputs: for every sample mentioned in "allClusterInfo.RDS", a file each for the DNA, RNA and protein levels of the contribution of each genome to the node is written to disk
#### "<sampleID>".<DNA/RNA/PROT>clusterPerNode.tsv"
#
#written by Anna Heintz-Buschart (July 2015)

nwNodes <- read.delim("150705_KOs_in_NW.tsv",stringsAsFactors=F)
clusterInfo <- readRDS("allClusterInfo.RDS")

library(rmongodb)

getExprData <- function(sampleID,cI=clusterInfo,nwAnno=nwNodes){
  #matches for pipelines:
  matchFun <- mongo.bson.from.JSON(paste('{"$match": {"genes.node": {"$exists": "true"}} }',sep=""))
  matchSample <- mongo.bson.from.JSON(paste('{"$match": {"sample": "',sampleID,'"} }',sep=""))
  matchProt <- mongo.bson.from.JSON(paste('{"$match": {"genes.proteinIdentification": "uniquely"} }',sep=""))
  #unwinds for pipelines:
  unwindGenes <- mongo.bson.from.JSON('{"$unwind": "$genes"}')
  unwindProteins <- mongo.bson.from.JSON('{"$unwind": "$genes.proteinIdentification"}')
  #grouping for pipelines:
  groupGeneA <- mongo.bson.from.JSON('{"$group": {"_id": "$_id", "gene": {"$push": "$genes.gene"},
                                     "aveCovDNA":{"$push":"$genes.aveCovDNA"},"aveCovRNA":{"$push":"$genes.aveCovRNAfw"},
                                     "fun":{"$push":"$genes.node"},"cluster":{"$push":"$cluster"}}}')
  #
  groupGeneProt <- mongo.bson.from.JSON('{"$group": {"_id": "$_id", "gene": {"$push": "$genes.gene"},
                                        "proteinIdentification":{"$push":"$genes.proteinIdentification"},
                                        "proteinArea":{"$push":"$genes.proteinArea"}}}')  
  #projections for pipelines:
  projectGeneA <- mongo.bson.from.JSON('{"$project": {"_id": 0, "fun":1, "gene":1,"cluster":1, "aveCovDNA":1,"aveCovRNA":1}}') 
  projectGeneProt <- mongo.bson.from.JSON('{"$project": {"_id": 0, "gene":1,"proteinIdentification":1,"proteinArea":1}}')
  
  if(!mongo.is.connected(mongo)) {
    stop("Mongo is not connected.")
  }else{
    #aggregation pipelines:
    genesA <- list()
    genesP <- list()
    for(clusOI in cI$cluster[cI$sample==sampleID&cI$cluster!="S"]){
      print(clusOI)
      matchCluster <- mongo.bson.from.JSON(paste('{"$match": {"cluster": "',clusOI,'"} }',sep=""))
      genesA <- append(genesA,mongo.bson.to.list(mongo.aggregation(mongo,coll,list(matchSample,matchCluster,matchFun,unwindGenes,matchFun,groupGeneA,projectGeneA)))$result)
      genesP <- append(genesP,mongo.bson.to.list(mongo.aggregation(mongo,coll,list(matchSample,matchCluster,matchFun,unwindGenes,matchFun,matchProt,groupGeneProt,projectGeneProt)))$result)
    }
    tFeat <- do.call(rbind,lapply(genesA,function(x) cbind(x$gene,x$fun,x$cluster,x$aveCovDNA,x$aveCovRNA)))
    tFeat[sapply(tFeat[,2],length)>1,2] <- sapply(tFeat[sapply(tFeat[,2],length)>1,2],function(x) paste(x,sep=";",collapse=";"))
    tFeat <- data.frame("gene"=unlist(tFeat[,1]),"node"=unlist(tFeat[,2]),"cluster"=unlist(tFeat[,3]),"aveCovDNA"=as.numeric(unlist(tFeat[,4])),
                        "aveCovRNA"=as.numeric(unlist(tFeat[,5])),stringsAsFactors=F)
    
    pFeat <- do.call(rbind,lapply(genesP,function(x) cbind(x$gene, x$proteinArea)))
    pFeat[sapply(pFeat[,2],length)>1,2] <- sapply(pFeat[sapply(pFeat[,2],length)>1,2], mean)
    pFeat <- data.frame("gene"=unlist(pFeat[,1]),"proteinArea"=as.numeric(unlist(pFeat[,2])),stringsAsFactors=F)
    aFeat <- merge(tFeat,pFeat,by="gene",all.x=T)
    aFeat[is.na(aFeat)] <- 0
    aFeat <- unique(aFeat)
    if(!"proteinArea" %in% colnames(aFeat)) aFeat$proteinArea <- 0
    
    dnaFeat <- tapply(aFeat$aveCovDNA,list(aFeat$cluster,aFeat$node),sum)
    dnaFeat[is.na(dnaFeat)] <- 0
    
    rnaFeat <- tapply(aFeat$aveCovRNA,list(aFeat$cluster,aFeat$node),sum)
    rnaFeat[is.na(rnaFeat)] <- 0
    
    protFeat <- tapply(aFeat$proteinArea,list(aFeat$cluster,aFeat$node),sum)
    protFeat[is.na(protFeat)] <- 0
    return(list(dnaFeat,rnaFeat,protFeat))
  } 
}

dir.create("./perSampleNodes2Pop/")
mongo <- mongo.create()
db <- "mydb" 
coll <- "mydb.must" 
for(sampleID in unique(clusterInfo$sample)){
  exTab <- getExprData(sampleID)
  
  dTab <- merge(exTab[[1]],clusterInfo[clusterInfo$sample==sampleID,c("uniqueEss","motuPresent","cluster")],by.y="cluster",by.x=0)
  rTab <- merge(exTab[[2]],clusterInfo[clusterInfo$sample==sampleID,c("uniqueEss","motuPresent","cluster")],by.y="cluster",by.x=0)
  pTab <- merge(exTab[[3]],clusterInfo[clusterInfo$sample==sampleID,c("uniqueEss","motuPresent","cluster")],by.y="cluster",by.x=0)
  
  dOut <- t(dTab[,-c(1,ncol(dTab)-0:1)][order(dTab$uniqueEss,decreasing=T),])
  colnames(dOut) <- dTab[,1][order(dTab$uniqueEss,decreasing=T)]
  rOut <- t(rTab[,-c(1,ncol(rTab)-0:1)][order(rTab$uniqueEss,decreasing=T),])
  colnames(rOut) <- rTab[,1][order(rTab$uniqueEss,decreasing=T)]
  pOut <- t(pTab[,-c(1,ncol(pTab)-0:1)][order(pTab$uniqueEss,decreasing=T),])
  colnames(pOut) <- pTab[,1][order(pTab$uniqueEss,decreasing=T)]
  
  write.table(dOut,paste("./perSampleNodes2Pop/",sampleID,".DNAclusterPerNode.tsv",sep=""),sep="\t",quote=F)
  write.table(rOut,paste("./perSampleNodes2Pop/",sampleID,".RNAclusterPerNode.tsv",sep=""),sep="\t",quote=F)
  write.table(pOut,paste("./perSampleNodes2Pop/",sampleID,".PROTclusterPerNode.tsv",sep=""),sep="\t",quote=F)
}
mongo.destroy(mongo)


