#this script runs to get all phylogenetic marker genes from mocat and rpoB from the MongoDB
#it writes a file with the gene names, samples and cluster IDs per gene class

# it should not be expected to run anywhere else other than the MuSt environment

#written by Anna Heintz-Buschart, August 2015

library(rmongodb)

getExprData <- function(marker,sampleID){
  #matches for pipelines:
  #matchGene <- mongo.bson.from.JSON(paste('{"$match": {"genes": {"$exists": "true"}} }',sep=""))
  matchSample <- mongo.bson.from.JSON(paste('{"$match": {"sample": "',sampleID,'"} }',sep=""))
  if(marker!="rpoB"){
    matchMarker <- mongo.bson.from.JSON(paste('{"$match": {"genes.mOTUbestMarkerGene": "',marker,'"} }',sep=""))
  }else{
    matchMarker <- mongo.bson.from.JSON(paste('{"$match": {"genes.essentialGene": "TIGR02013"} }',sep=""))
  }
  matchComp <- mongo.bson.from.JSON(paste('{"$match": {"genes.completeness": "complete"} }',sep=""))
  #unwinds for pipelines:
  unwindGenes <- mongo.bson.from.JSON('{"$unwind": "$genes"}')
  #grouping for pipelines:
  if(marker!="rpoB"){
    groupGeneA <- mongo.bson.from.JSON('{"$group": {"_id": "$_id", "gene": {"$push": "$genes.gene"},
                                     "fun":{"$push":"$genes.mOTUbestMarkerGene"},"sample":{"$push":"$sample"},
                                       "cluster":{"$push":"$cluster"}, "completeness":{"$push":"$genes.completeness"}}}')
  }else{
    groupGeneA <- mongo.bson.from.JSON('{"$group": {"_id": "$_id", "gene": {"$push": "$genes.gene"},
                                     "fun":{"$push":"$genes.essentialGene"},"sample":{"$push":"$sample"},
                                       "cluster":{"$push":"$cluster"}, "completeness":{"$push":"$genes.completeness"}}}')
  }
  #projections for pipelines:
  projectGeneA <- mongo.bson.from.JSON('{"$project": {"_id": 0, "fun":1, "gene":1,"cluster":1, "sample":1,"completeness":1}}') 
  
  if(!mongo.is.connected(mongo)) {
    stop("Mongo is not connected.")
  }else{
    #aggregation pipelines:
    genesA <- mongo.bson.to.list(mongo.aggregation(mongo,coll,
                                                   list(matchSample,matchMarker,unwindGenes,matchMarker,matchComp,groupGeneA,projectGeneA)))$result

    aFeat <- do.call(rbind,lapply(genesA,function(x) cbind(x$gene,x$fun,x$cluster,x$sample,x$completeness)))
    aFeat <- data.frame("gene"=unlist(aFeat[,1]),"anno"=unlist(aFeat[,2]),"cluster"=unlist(aFeat[,3]),"sample"=unlist(aFeat[,4]),
                        "completeness"=unlist(aFeat[,5]),stringsAsFactors=F)
    aFeat <- unique(aFeat)
    return(aFeat[,c("gene","sample","cluster")])
  } 
}

clusterInfo <- readRDS("../Bmaps/allClusterInfo.RDS")
allS <- sort(unique(clusterInfo$sample))
ids <- read.delim("combinedIds",header=F,stringsAsFactors=F)
colnames(ids) <- c("fam","sample")

mongo <- mongo.create()
db <- "mydb" 
coll <- "mydb.must" 
for(sample in allS){
  fam <- ids$fam[ids$sample==sample]
  for(mark in c("rpoB","COG0012","COG0016","COG0018","COG0172","COG0215","COG0495","COG0525","COG0533","COG0541","COG0552")){
    print(paste(sample,mark))
    feats <- getExprData(mark,sample)
    write.table(feats,paste("/work/projects/ecosystem_biology/MUST/CombinedAssembly/Annotation/mOTUgenes/",fam,"/",sample,"/",mark,".","geneNamesClusters.tsv",sep=""),
                col.names=F,row.names=F,quote=F,sep="\t")
  }
}
mongo.destroy(mongo)
