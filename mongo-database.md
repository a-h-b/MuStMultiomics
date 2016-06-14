For the faecal samples of the MuSt that were analyzed on the metagenomic, metatranscriptomic and metaproteomic levels, we assembled a total of more than 26 million contigs containing approximately 30 million genes. For each of these contigs, we had interesting data, such as length, metagenomic coverage, position and coverage of variants, membership in a bin, position on a BH-SNE map, putative taxonomy and of course the genes they harbour. For the genes, we had similar data, such as the position and sense on the contig, coverage by metagenomic and metatranscriptomic reads, detection of the protein product, functional annotation(s), essentiality, and taxonomic annotations. 

My favourite way of visualizing and statistically analyzing data is by using [R](https://www.r-project.org/). But even sample-wise datasets are huge for R to store, search and filter. Of course, there are some summarizing data sets, such as the abundance of every mOTU in every sample (matrix of 500 x 36) or the expression of every KO (matrix of a few thousand x 36) in every sample which can be easily worked with in R. But to find for example the taxonomy of every gene annotated as K00001 in every sample, along with its expression on the metatranscriptomic and metaproteomic level would be very difficult (or better: time consuming) using pure R. In addition, someone else might want to use a different tool to analyze the data.

I finally decided to build a MongoDB database (https://www.mongodb.org/), mainly because the contigs and genes naturally form a nested structure that is well representable by Mongo documents and because I hope that Mongo would scale well for adaptation to a larger project with a biologically meaningful sample size. I was also intrigued by the fact that it has a 2D indexing which can be used for searching, but I have actually never gotten around to implementing using the BH-SNE maps as interface to the database. MongoDBs can be accessed with a large number of tools or programming languages (I have only used R and python here) and it is nice and fast.

A python script (final version is [`150928_mongifyMust.py`](150928_mongifyMust.py)) was used to fill the database with most of the metaG, metaT and metaP data that was created within the MuSt with the exception of the actual sequences. This script is made for the MuSt with its particular system of sorting data by families and then samples. It is made up mostly of paths to the files that were created one by one without the use of the database. Therefore reusing this script in any other environment is not recommended. The resulting structure with the names of the fields is displayed here for your reference if you want to use some of the code to access it from R:
[database structure](http://git-r3lab.uni.lu/anna.buschart/MuStMultiomics/blob/master/figS11_DB.png)

The script [`151020_funOIMongoWS.R`](151020_funOIMongoWS.R) exemplifies how the database can be accessed from R to retrieve expression data for genes with a function of interest and the contig bins these genes belong to. To just give the most important points, here is what always needs to be done to interact with a MongoDB using the rmongodb-package (https://cran.r-project.org/web/packages/rmongodb/index.html).

```
library(rmongodb)
mongo <- mongo.create()
db <- "mydb" 
coll <- "mydb.must" 
# access the database here using the custom function getExprTab and the function of interest:
exTab <- getExprTab(funOI)
mongo.destroy(mongo)
```

The function `getExprTab()` then fetches the data into two lists via two aggregation pipelines:

```
unwindGenes <- mongo.bson.from.JSON('{"$unwind": "$genes"}')
matchAnno <- mongo.bson.from.JSON(paste('{"$match": {"genes.bestAnnotation": "',funOI,'"} }',sep=""))
matchProt <- mongo.bson.from.JSON(paste('{"$match": {"genes.proteinIdentification": "uniquely"} }',sep=""))
groupDNARNAname <- mongo.bson.from.JSON('{"$group": {"_id": "$_id", "aveCovDNA": {"$push": "$genes.aveCovDNA"},"sample":{"$push":"$sample"},"gene":{"$push":"$genes.gene"},"cluster":{"$push": "$cluster"},"aveCovRNA":{"$push": "$genes.aveCovRNAfw"}}}')
projectDNARNAname <- mongo.bson.from.JSON('{"$project": {"_id": 0, "gene":1,"aveCovDNA": 1, "cluster" :1,"aveCovRNA":1,"sample":1}}')
groupProtname <- mongo.bson.from.JSON('{"$group": {"_id": "$_id", "gene": {"$push": "$genes.gene"},"sample":{"$push":"$sample"},"proteinIdentification":{"$push":"$genes.proteinIdentification"},"proteinArea":{"$push":"$genes.proteinArea"}}}') 
projectProtname <- mongo.bson.from.JSON('{"$project": {"_id": 0, "gene":1,"proteinIdentification":1,"proteinArea":1,"sample":1}}')
genesA <- mongo.bson.to.list(mongo.aggregation(mongo,coll,list(matchAnno,unwindGenes,matchAnno,groupDNARNAname,projectDNARNAname)))$result
genesP <- mongo.bson.to.list(mongo.aggregation(mongo,coll,list(matchAnno,unwindGenes,matchAnno,matchProt,groupProtname,projectProtname)))$result
```
As not all genes are found on the protein level, the proteins have their own pipeline.
The list with the data that is present for all genes is transformed into a data.frame like this:

```
res <- do.call(rbind,lapply(genesA,function(x) cbind(x$sample,x$gene,x$cluster,x$aveCovDNA,x$aveCovRNA)))
res <- data.frame("sample"=res[,1],"gene"=res[,2],"cluster"=res[,3],"aveCovDNA"=as.numeric(res[,4]),"aveCovRNA"=as.numeric(res[,5]),stringsAsFactors=F)
```
If there is data for the protein level, this is merged to the data.frame, otherwise 0 is returned for the protein abundance:

```
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
```

Another example is given in [`150928_MUST_relatedClusterWSFromMongo.R`](150928_MUST_relatedClusterWSFromMongo.R). Here, a function is used which returns all genes of one bin ("cluster", clusOI) of contigs in a sample (sampleID) with their functions and levels on the different omic levels. This again uses aggregation pipelines with the following steps:

```
matchFun <- mongo.bson.from.JSON(paste('{"$match": {"genes.bestAnnotation": {"$exists": "true"}} }',sep=""))
matchCluster <- mongo.bson.from.JSON(paste('{"$match": {"cluster": "',clusOI,'"} }',sep=""))
matchSample <- mongo.bson.from.JSON(paste('{"$match": {"sample": "',sampleID,'"} }',sep=""))
matchProt <- mongo.bson.from.JSON(paste('{"$match": {"genes.proteinIdentification": "uniquely"} }',sep=""))
unwindGenes <- mongo.bson.from.JSON('{"$unwind": "$genes"}')
unwindProteins <- mongo.bson.from.JSON('{"$unwind": "$genes.proteinIdentification"}')
groupGeneA <- mongo.bson.from.JSON('{"$group": {"_id": "$_id", "gene": {"$push": "$genes.gene"}, "aveCovDNA":{"$push":"$genes.aveCovDNA"}, "aveCovRNA":{"$push":"$genes.aveCovRNAfw"}, "readsRNA":{"$push":"$genes.readsRNAfw"},"fun":{"$push":"$genes.bestAnnotation"}}}')
groupGeneProt <- mongo.bson.from.JSON('{"$group": {"_id": "$_id", "gene": {"$push": "$genes.gene"}, "proteinIdentification":{"$push":"$genes.proteinIdentification"}, "proteinArea":{"$push":"$genes.proteinArea"}}}')  
projectGeneA <- mongo.bson.from.JSON('{"$project": {"_id": 0, "fun":1, "gene":1, "aveCovDNA":1,"aveCovRNA":1,"readsRNA":1}}')
projectGeneProt <- mongo.bson.from.JSON('{"$project": {"_id": 0, "gene":1,"proteinIdentification":1,"proteinArea":1}}')
```
The results are fed into lists.

```
genesA <- mongo.bson.to.list(mongo.aggregation(mongo,coll,list(matchCluster,matchSample,unwindGenes,matchFun,groupGeneA,projectGeneA)))$result
genesP <- mongo.bson.to.list(mongo.aggregation(mongo,coll,list(matchCluster,matchSample,unwindGenes,matchFun,matchProt,groupGeneProt,projectGeneProt)))$result
```
Which are again converted to data.frames using `do.call`, `rbind` and `lapply`:

```
tFeat <- do.call(rbind,lapply(genesA,function(x) cbind(x$gene,x$fun,x$aveCovDNA,x$aveCovRNA,x$readsRNA)))
tFeat[sapply(tFeat[,2],length)>1,2] <- sapply(tFeat[sapply(tFeat[,2],length)>1,2],function(x) paste(x,sep=";",collapse=";"))
tFeat <- data.frame("gene"=unlist(tFeat[,1]),"bestAnno"=unlist(tFeat[,2]),"aveCovDNA"=as.numeric(unlist(tFeat[,3])),"aveCovRNA"=as.numeric(unlist(tFeat[,4])),"readsRNA"=as.numeric(unlist(tFeat[,5])),stringsAsFactors=F)
```

The scripts [`virusGenesMongo.R`](virusGenesMongo.R) and [`eukaryoticGenesMongo.R`](eukaryoticGenesMongo.R) exemplify how functional annotations and coverage data for genes can be accessed by the gene ID and the sample ID.

Similarly, the contribution of all binned genome-level populations to the abundance of nodes in a reconstructed network can be retrieved using the script [`getNWexprMongoAllSamples.R`](getNWexprMongoAllSamples.R).


