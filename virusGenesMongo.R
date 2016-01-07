# this script contains code used in the MuSt project to find functional annotations in the MongoDB, in this case for viral genes 

#inputs: uses - "MGRastIDs.txt", a file linking MG-RAST IDs to the IDs used throughout the MuSt
#             - for each sample, the output of "MGRASTgeneLevelTax.R" for the viral genes  
#output: an R workspace containing most importantly two tables "allVirFunc" and "allVir", which contain functional and 
### coverage data and taxonomic annotations from Kraken for all viral genes in all samples

# this script is summarized from several scripts that were used in the MuSt to get the same data
#written by Anna Heintz-Buschart (Nov-Dec 2015, summarized Jan 2016)

library(rmongodb)

MGRASTids <- read.delim("MGRastIDs.txt",header=F,stringsAsFactors=F)
colnames(MGRASTids) <- c("fam","sample","mgrid")

geneFiles <- vector()
for(FAM in unique(MGRASTids$fam)){
  geneFiles <- append(geneFiles,paste("../MGRAST/",FAM,"/",list.files(path=paste("../MGRAST/",FAM,sep=""),pattern="viruses.tsv"),sep=""))
}
first <- T
for(gf in geneFiles){
  LIB <- gsub("^../MGRAST/must_m_0./","",gsub(".vir.+","",gf))
  if(first){
    genesTab <- data.frame("sample"=LIB,read.delim(gf,stringsAsFactors=F),stringsAsFactors=F)
    first <- F
  }else{
    genesTab <- rbind(genesTab,data.frame("sample"=LIB,read.delim(gf,stringsAsFactors=F),stringsAsFactors=F))
  }
}


mongo <- mongo.create()
db <- "mydb" 
coll <- "mydb.must" 

if(mongo.is.connected(mongo)) {
  matchProt <- mongo.bson.from.JSON(paste('{"$match": {"genes.proteinIdentification": "uniquely"} }',sep=""))
  unwindGenes <- mongo.bson.from.JSON('{"$unwind": "$genes"}')
  groupDNARNAname <- mongo.bson.from.JSON('{"$group": {"_id": "$_id", "aveCovDNA": {"$push": "$aveCov"},"sample":{"$push":"$sample"},
                                          "gene":{"$push":"$genes.gene"},"cluster":{"$push": "$cluster"},"aveCovRNA":{"$push": "$genes.aveCovRNAfw"},
                                          "readCovRNA":{"$push":"$genes.readsRNAfw"},"krakenSpecies":{"$push":"$krakenSpecies"},
                                          "krakenFamily":{"$push":"$krakenFamily"},"krakenOrder":{"$push":"$krakenOrder"},"krakenKingdom":{"$push":"$krakenKingdom"}}}')
  groupFunc <- mongo.bson.from.JSON('{"$group": {"_id": "$_id", "bestAnno": {"$push": "$genes.bestAnnotation"},"sample":{"$push":"$sample"},
                                    "gene":{"$push":"$genes.gene"},"length":{"$push": "$length"}}}')
  groupProtname <- mongo.bson.from.JSON('{"$group": {"_id": "$_id", "gene": {"$push": "$genes.gene"},"sample": {"$push": "$sample"},"proteinIdentification":{"$push":"$genes.proteinIdentification"},"proteinArea":{"$push":"$genes.proteinArea"}}}') 
  projectDNARNAname <- mongo.bson.from.JSON('{"$project": {"_id": 0, "gene":1,"aveCovDNA": 1, "cluster" :1,"aveCovRNA":1,"readCovRNA":1,"sample":1,"krakenSpecies":1,"krakenFamily":1,"krakenOrder":1,"krakenKingdom":1}}')
  projectFunc <- mongo.bson.from.JSON('{"$project": {"_id": 0, "gene":1,"bestAnno": 1, "sample":1,"length":1}}')
  projectProtname <- mongo.bson.from.JSON('{"$project": {"_id": 0, "gene":1,"proteinIdentification":1,"proteinArea":1,"sample":1}}')
  
  genesA <- list()
  genesF <- list()
  genesP <- list()
  for(i in 1:nrow(genesTab)){
    matchGene <- mongo.bson.from.JSON(paste('{"$match": {"genes.gene": "',genesTab$gene[i],'"} }',sep=""))
    matchSample <- mongo.bson.from.JSON(paste('{"$match": {"sample": "',genesTab$sample[i],'"} }',sep=""))
    genesF <- append(genesF,mongo.bson.to.list(mongo.aggregation(mongo,coll,list(matchGene,matchSample,unwindGenes,matchGene,groupFunc,projectFunc)))$result)
    genesA <- append(genesA,mongo.bson.to.list(mongo.aggregation(mongo,coll,list(matchGene,matchSample,unwindGenes,matchGene,groupDNARNAname,projectDNARNAname)))$result)
    genesP <- append(genesP,mongo.bson.to.list(mongo.aggregation(mongo,coll,list(matchGene,matchSample,unwindGenes,matchGene,matchProt,groupProtname,projectProtname)))$result)
    if(length(genesA)==0){ 
      warning("No genes found")
      return(NULL)
    }
  }
  }else {
    stop("Mongo is not connected.")
  }
mongo.destroy(mongo)
system("bash stopMongo.sh")

fFeat <- do.call(rbind,lapply(genesF,function(x) {
  if(length(x$bestAnno)>0) bA <- paste(x$bestAnno,sep=";",collapse=";") else bA <- ""
  cbind(x$gene,x$sample,bA,x$length)}))

fFeat <- data.frame("gene"=unlist(fFeat[,1]),"sample"=unlist(fFeat[,2]),"bestAnno"=unlist(fFeat[,3]),
                    "contiglength"=as.numeric(unlist(fFeat[,4])),stringsAsFactors=F)    

fFeat[is.na(fFeat)] <- 0
fFeat <- unique(fFeat)
allVirFunc <- merge(genesTab,fFeat,by.x=c(1,2),by.y=c(2,1),all.x=T)

####################

tFeat <- do.call(rbind,lapply(genesA,function(x){
  if(length(x$krakenSpecies)==0) ks <- "unknown" else ks <- x$krakenSpecies
  if(length(x$krakenFamily)==0) kf <- "unknown" else kf <- x$krakenFamily
  if(length(x$krakenOrder)==0) ko <- "unknown" else ko <- x$krakenOrder
  if(length(x$krakenKingdom)==0) kk <- "unknown" else kk <- x$krakenKingdom
  cbind(x$gene,x$sample,x$aveCovDNA,x$aveCovRNA,x$readCovRNA,x$cluster,ks,kf,ko,kk)
}))

tFeat <- data.frame("gene"=unlist(tFeat[,1]),"sample"=unlist(tFeat[,2]),"aveCovDNA"=as.numeric(unlist(tFeat[,3])),
                    "aveCovRNA"=as.numeric(unlist(tFeat[,4])),"RNAreads"=as.numeric(unlist(tFeat[,5])),"cluster"=unlist(tFeat[,6]),
                    "krakenSpecies"=unlist(tFeat[,7]),"krakenFamily"=unlist(tFeat[,8]),"krakenOrder"=unlist(tFeat[,9]),"krakenKingdom"=unlist(tFeat[,10]),stringsAsFactors=F)    

pFeat <- do.call(rbind,lapply(genesP,function(x) cbind(x$gene, x$sample, x$proteinArea)))
pFeat[sapply(pFeat[,3],length)>1,3] <- sapply(pFeat[sapply(pFeat[,3],length)>1,3], mean)
pFeat <- data.frame("gene"=unlist(pFeat[,1]),"sample"=unlist(pFeat[,2]),"proteinArea"=as.numeric(unlist(pFeat[,3])),stringsAsFactors=F)
if(nrow(pFeat)>0) aFeat <- merge(tFeat,pFeat,by=c(1,2),all.x=T) else aFeat <- tFeat
aFeat[is.na(aFeat)] <- 0
aFeat <- unique(aFeat)
if(!"proteinArea" %in% colnames(aFeat)) aFeat$proteinArea <- 0   
allVir <- merge(genesTab,aFeat,by.x=c(1,2),by.y=c(2,1),all.x=T)


save.image("MGRASTvirFuncsCovWS.Rdata")

