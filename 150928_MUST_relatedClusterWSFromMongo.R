#this script gets all genes of related genomes, along with their expression levels on the different omic levels and the best annotation
#### reads are summed up by functions and DESeq-normalized
#### and creates a number of workspace with the data

# the script uses allClusterInfo.RDS, which contains statistics on all bins of all samples
# the function getExprData() defined in lines 23-68 does the interfacing with MongoDB
# the loop in line 76-109 does the gathering by same taxonomic annotation

# written by Anna Heintz-Buschart, September 2015 - this version is from January 2016, shortened to the part that was actually used in the MuSt

clusterInfo <- readRDS("allClusterInfo.RDS")

library(rmongodb)
mongo <- mongo.create()
db <- "mydb"
coll <- "mydb.must"
if(mongo.is.connected(mongo)) print(paste("found connection to mongod -",mongo.count(mongo,coll),"documents in collection"))

library(DESeq2)

#####

getExprData <- function(sampleID,clusOI){
  #matches for pipelines:
  matchFun <- mongo.bson.from.JSON(paste('{"$match": {"genes.bestAnnotation": {"$exists": "true"}} }',sep=""))
  matchCluster <- mongo.bson.from.JSON(paste('{"$match": {"cluster": "',clusOI,'"} }',sep=""))
  matchSample <- mongo.bson.from.JSON(paste('{"$match": {"sample": "',sampleID,'"} }',sep=""))
  matchProt <- mongo.bson.from.JSON(paste('{"$match": {"genes.proteinIdentification": "uniquely"} }',sep=""))
  #unwinds for pipelines:
  unwindGenes <- mongo.bson.from.JSON('{"$unwind": "$genes"}')
  unwindProteins <- mongo.bson.from.JSON('{"$unwind": "$genes.proteinIdentification"}')
  #grouping for pipelines:
  groupGeneA <- mongo.bson.from.JSON('{"$group": {"_id": "$_id", "gene": {"$push": "$genes.gene"},
                                     "aveCovDNA":{"$push":"$genes.aveCovDNA"},"aveCovRNA":{"$push":"$genes.aveCovRNAfw"},
                                     "readsRNA":{"$push":"$genes.readsRNAfw"},"fun":{"$push":"$genes.bestAnnotation"}}}')
  groupGeneProt <- mongo.bson.from.JSON('{"$group": {"_id": "$_id", "gene": {"$push": "$genes.gene"},
                                        "proteinIdentification":{"$push":"$genes.proteinIdentification"},
                                        "proteinArea":{"$push":"$genes.proteinArea"}}}')  
  #projections for pipelines:
  projectGeneA <- mongo.bson.from.JSON('{"$project": {"_id": 0, "fun":1, "gene":1, "aveCovDNA":1,"aveCovRNA":1,"readsRNA":1}}')
  projectGeneProt <- mongo.bson.from.JSON('{"$project": {"_id": 0, "gene":1,"proteinIdentification":1,"proteinArea":1}}')
  
  if(!mongo.is.connected(mongo)) {
    stop("Mongo is not connected.")
  }else{
    #aggregation pipelines:
    genesA <- mongo.bson.to.list(mongo.aggregation(mongo,coll,list(matchCluster,matchSample,unwindGenes,matchFun,groupGeneA,
                                                                   projectGeneA)))$result
    if(length(genesA)==0){ 
      warning("No genes found")
    } else {
      tFeat <- do.call(rbind,lapply(genesA,function(x) cbind(x$gene,x$fun,x$aveCovDNA,x$aveCovRNA,x$readsRNA)))
      tFeat[sapply(tFeat[,2],length)>1,2] <- sapply(tFeat[sapply(tFeat[,2],length)>1,2],function(x) paste(x,sep=";",collapse=";"))
      tFeat <- data.frame("gene"=unlist(tFeat[,1]),"bestAnno"=unlist(tFeat[,2]),"aveCovDNA"=as.numeric(unlist(tFeat[,3])),
                          "aveCovRNA"=as.numeric(unlist(tFeat[,4])),"readsRNA"=as.numeric(unlist(tFeat[,5])),stringsAsFactors=F)
      
      genesP <- mongo.bson.to.list(mongo.aggregation(mongo,coll,list(matchCluster,matchSample,unwindGenes,matchFun,matchProt,groupGeneProt,
                                                                     projectGeneProt)))$result
      pFeat <- do.call(rbind,lapply(genesP,function(x) cbind(x$gene, x$proteinArea)))
      pFeat[sapply(pFeat[,2],length)>1,2] <- sapply(pFeat[sapply(pFeat[,2],length)>1,2], mean)
      pFeat <- data.frame("gene"=unlist(pFeat[,1]),"proteinArea"=as.numeric(unlist(pFeat[,2])),stringsAsFactors=F)
      aFeat <- merge(tFeat,pFeat,by=1,all.x=T)
      aFeat[is.na(aFeat)] <- 0
      if(!"proteinArea" %in% colnames(aFeat)) aFeat$proteinArea <- 0
      return(aFeat)
    } 
  }
}


motuDup <- names(table(sapply(grep(";",clusterInfo$motuPresent,invert=T,value=T),
                              function(x)unlist(strsplit(x,"\\("))[1]))
                 [table(sapply(grep(";",clusterInfo$motuPresent,invert=T,value=T),
                               function(x)unlist(strsplit(x,"\\("))[1]))>1])

for(currMotu in motuDup){
  clOI <- clusterInfo[grepl(paste(currMotu,"\\(",sep=""),clusterInfo$motuPresent)&!grepl(";",clusterInfo$motuPresent),c("cluster","sample")]
  first <- T
  for(i in 1:nrow(clOI)){
    sI <- clOI$sample[i]
    clI <- clOI$cluster[i]
    if(!is.na(sI)&!is.na(clI)){
      cwd <- paste(sI,clI,sep="_")
      if(first){
        exprall <- data.frame("cluster"=cwd,getExprData(sI,clI),stringsAsFactors=F)
        first <- F
      }else{
        texp <- data.frame("cluster"=cwd,getExprData(sI,clI),stringsAsFactors=F)
        exprall <- rbind(exprall,texp)
      }
    }
  }
  pres <- tapply(exprall$gene,list(exprall$bestAnno,exprall$cluster),function(x)as.numeric(length(x)>0))
  pres[is.na(pres)] <-0
  if(ncol(pres)>1&nrow(pres)>1){    
    readR <- tapply(exprall$readsRNA,list(exprall$bestAnno,exprall$cluster),sum)
    readR[is.na(readR)] <-0

    gps <- c(grep("P",colnames(pres)),grep("G",colnames(pres)))
    if(length(gps)>2){
      readRA <- readR[,gps]
      readRA <- readRA[apply(dcov[,gps],1,function(x)all(x>0)),]
      readDDS <- DESeqDataSetFromMatrix(readRA,data.frame("sample"=colnames(readRA)),~sample)
      ts2 <- estimateSizeFactors(readDDS)
      readRN <- counts(ts2,normalized=T)
      save.image(paste(currMotu,"_WS.Rdata",sep=""))
    } 
  }
}

mongo.destroy(mongo)
