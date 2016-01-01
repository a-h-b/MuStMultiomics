#functions
extrCont <- function(geneNames) sub("gene","",gsub("_r?[[:digit:]]+$","",geneNames))

#arguments from command line
args<-commandArgs(TRUE)
FAM <- args[1]
LIB <- args[2]

#translation of IDs
ids <- read.delim(paste("/work/projects/ecosystem_biology/MUST/CombinedAssembly/Annotation/sequenceInfo/",FAM,"/ids",sep=""),stringsAsFactors=F,header=F)
visits <- read.delim(paste("/work/projects/ecosystem_biology/MUST/CombinedAssembly/Annotation/mapping/",FAM,"/visits",sep=""),stringsAsFactors=F,header=F)

#genewise table (all genes)
###
protInf <- read.delim(paste("/work/projects/ecosystem_biology/MUST/CombinedAssembly/Annotation/sequenceInfo/",FAM,"/",LIB,"/gene.prediction.assembly.Prodigal.500/contig.Prodigal.tab",sep=""),stringsAsFactors=F,header=F)
colnames(protInf) <- c("gene","sense","start","end","length","startCodon","stopCodon","completeness")
protInf$contig <- extrCont(protInf$gene)
protInf$kind <- rep("protein",nrow(protInf))
rnaInf <- read.delim(paste("/work/projects/ecosystem_biology/MUST/CombinedAssembly/Annotation/hmms/",FAM,"/",LIB,"/rRNAgenes.bac.tab",sep=""),stringsAsFactors=F)
rnaInf$startCodon <- rep("NA",nrow(rnaInf))
rnaInf$stopCodon <- rep("NA",nrow(rnaInf))
geneInf <- rbind(protInf[,c("gene","contig","sense","start","end","length","startCodon","stopCodon","completeness","kind")],rnaInf[,c("gene","contig","sense","start","end","length","startCodon","stopCodon","completeness","kind")])
# DNA coverage of genes:
dLIB <- ids$V1[ids$V2==LIB]
cov <- read.delim(paste("/work/projects/ecosystem_biology/MUST/CombinedAssembly/Annotation/mapping/",FAM,"/",LIB,"/",dLIB,".DNAonGenesrRNA.cov.tsv",sep=""),stringsAsFactors=F,colClasses=c("character","NULL","numeric","NULL"))
geneInf <- merge(geneInf,cov,by.x=1,by.y=1,all.x=T)
colnames(geneInf)[ncol(geneInf)] <- "aveCovDNA"
geneInf$aveCovDNA[is.na(geneInf$aveCovDNA)] <- 0
# RNA coverage of genes:
rnafile <- paste("/work/projects/ecosystem_biology/MUST/CombinedAssembly/Annotation/mapping/",FAM,"/",LIB,"/",LIB,".RNAonGenesrRNA.cov.tsv",sep="")
cov <- read.delim(rnafile,stringsAsFactors=F,colClasses=c("character","NULL","numeric","NULL"))
geneInf <- merge(geneInf,cov,by.x=1,by.y=1,all.x=T)
colnames(geneInf)[ncol(geneInf)] <- "aveCovRNAfw"
geneInf$aveCovRNAfw[is.na(geneInf$aveCovRNAfw)] <- 0
rnafile <- paste("/work/projects/ecosystem_biology/MUST/CombinedAssembly/Annotation/mapping/",FAM,"/",LIB,"/",LIB,".revcompRNAonGenesrRNA.cov.tsv",sep="")
cov <- read.delim(rnafile,stringsAsFactors=F,colClasses=c("character","NULL","numeric","NULL"))
geneInf <- merge(geneInf,cov,by.x=1,by.y=1,all.x=T)
colnames(geneInf)[ncol(geneInf)] <- "aveCovRNArc"
geneInf$aveCovRNArc[is.na(geneInf$aveCovRNArc)] <- 0
#variants on genes
var <- read.table(paste("/work/projects/ecosystem_biology/MUST/CombinedAssembly/Annotation/variants/",FAM,"/",LIB,"/variantsAnnotation.tsv",sep=""),header=T,stringsAsFactors=F)
gvar <- aggregate(var$pos,list(var$gene),length) 
geneInf <- merge(geneInf,gvar,by.x="gene",by.y=1,all.x=T)
colnames(geneInf)[ncol(geneInf)] <- "varPerMB"
geneInf$varPerMB[is.na(geneInf$varPerMB)] <- 0
geneInf$varPerMB <- geneInf$varPerMB*1000000/geneInf$length


#essential genes
essGenes <- read.table(paste("/work/projects/ecosystem_biology/MUST/CombinedAssembly/Annotation/hmms/",FAM,"/",LIB,"/contig.Prodigal.hmm.orfs.hits",sep=""),header=F,stringsAsFactors=F,quote="",strip.white=T,skip=3,colClasses=c("character","NULL","character","NULL","numeric",rep("NULL",14)))
dupliEss <- unique(names(table(essGenes$V1)[table(essGenes$V1)>1]))
dupliID <- dupliEss
for(dupGenes in dupliEss) dupliID[which(dupliID==dupGenes)] <- essGenes$V3[essGenes$V1==dupGenes][which.min(essGenes$V5[essGenes$V1==dupGenes])]
essGenes <- rbind(essGenes[!(essGenes$V1 %in% dupliEss),-3],data.frame("V1"=dupliEss,"V3"=dupliID,stringsAsFactors=F))
geneInf <- merge(geneInf,essGenes,by.x=1,by.y=1,all.x=T)
colnames(geneInf)[ncol(geneInf)] <- "essentialGene"
geneInf$essentialGene[is.na(geneInf$essentialGene)] <- "notEssential"
#kegg annotation
kegg <- read.delim(paste("/work/projects/ecosystem_biology/MUST/CombinedAssembly/Annotation/hmms/",FAM,"/",LIB,"/contig.Prodigal.faakegg.hits.tsv",sep=""),stringsAsFactors=F,colClasses=rep(c("character","NULL"),each=2))
geneInf <- merge(geneInf,kegg,by.x=1,by.y=1,all.x=T)
geneInf$KO[is.na(geneInf$KO)] <- "unknown"
#MetaCyc annotation
metCy <- read.delim(paste("/work/projects/ecosystem_biology/MUST/CombinedAssembly/Annotation/hmms/",FAM,"/",LIB,"/contig.Prodigal.faametacyc.hits.tsv",sep=""),stringsAsFactors=F,colClasses=rep(c("character","NULL"),each=2))
geneInf <- merge(geneInf,metCy,by.x=1,by.y=1,all.x=T)
geneInf$metaCycID[is.na(geneInf$metaCycID)] <- "unknown"
#swissprot annotation
sw <- read.delim(paste("/work/projects/ecosystem_biology/MUST/CombinedAssembly/Annotation/hmms/",FAM,"/",LIB,"/contig.Prodigal.faaswissprot.hits.tsv",sep=""),stringsAsFactors=F,colClasses=rep(c("character","NULL"),each=2))
geneInf <- merge(geneInf,sw,by.x=1,by.y=1,all.x=T)
geneInf$swissprotEC[is.na(geneInf$swissprotEC)] <- "unknown"
#pfam annotation
pfam <- read.delim(paste("/work/projects/ecosystem_biology/MUST/CombinedAssembly/Annotation/hmms/",FAM,"/",LIB,"/contig.Prodigal.faapfam.hits.tsv",sep=""),stringsAsFactors=F,colClasses=rep(c("character","NULL"),each=2))
geneInf <- merge(geneInf,pfam,by.x=1,by.y=1,all.x=T)
geneInf$pfamID[is.na(geneInf$pfamID)] <- "unknown"
#tigrPfam annotation
tigr <- read.delim(paste("/work/projects/ecosystem_biology/MUST/CombinedAssembly/Annotation/hmms/",FAM,"/",LIB,"/contig.Prodigal.faatigrpfam.hits.tsv",sep=""),stringsAsFactors=F,colClasses=rep(c("character","NULL"),each=2))
geneInf <- merge(geneInf,tigr,by.x=1,by.y=1,all.x=T)
geneInf$tigrID[is.na(geneInf$tigrID)] <- "unknown"

#amphora annotation
amp <- read.delim(paste("/work/projects/ecosystem_biology/MUST/CombinedAssembly/Annotation/hmms/",FAM,"/",LIB,"/amphora2.corrected.tsv",sep=""),stringsAsFactors=F)
geneInf <- merge(geneInf,amp,by.x=1,by.y=1,all.x=T)
for(i in ncol(geneInf)-0:7){
	colnames(geneInf)[i] <- paste(colnames(geneInf)[i],"Amphora",sep=".")
	geneInf[is.na(geneInf[,i]),i] <- "unknown"
	geneInf[geneInf[,i]=="",i] <- "unknown"
}
#mOTU genes (best hit)
mphFiles <- list.files(paste("/work/projects/ecosystem_biology/MUST/CombinedAssembly/Annotation/mOTUgenes/",FAM,"/",LIB,"/",sep=""),"bestHitPhylogeny")
mge <- read.delim(paste("/work/projects/ecosystem_biology/MUST/CombinedAssembly/Annotation/mOTUgenes/",FAM,"/",LIB,"/",mphFiles[1],sep=""),header=T,stringsAsFactors=F)
for(m in mphFiles[-1]){
mge <- rbind(mge,read.delim(paste("/work/projects/ecosystem_biology/MUST/CombinedAssembly/Annotation/mOTUgenes/",FAM,"/",LIB,"/",m,sep=""),header=T,stringsAsFactors=F))
}
geneInf <- merge(geneInf,mge[,c(4,5,15,17:23,1)],by.x="gene",by.y="queryID",all.x=T)
for(i in ncol(geneInf)-0:9){
        colnames(geneInf)[i] <- paste(colnames(geneInf)[i],"mOTUbest",sep=".")
        geneInf[is.na(geneInf[,i]),i] <- "unknown"
}
#mOTU genes (best hit out of present mOUTs (presence defined by read level analyis))
mphFiles <- list.files(paste("/work/projects/ecosystem_biology/MUST/CombinedAssembly/Annotation/mOTUgenes/",FAM,"/",LIB,"/",sep=""),"bestPresentHitPhylogeny")
mge <- read.delim(paste("/work/projects/ecosystem_biology/MUST/CombinedAssembly/Annotation/mOTUgenes/",FAM,"/",LIB,"/",mphFiles[1],sep=""),header=T,stringsAsFactors=F)
for(m in mphFiles[-1]){
mge <- rbind(mge,read.delim(paste("/work/projects/ecosystem_biology/MUST/CombinedAssembly/Annotation/mOTUgenes/",FAM,"/",LIB,"/",m,sep=""),header=T,stringsAsFactors=F))
}
geneInf <- merge(geneInf,mge[,c(4,5,15,17:23,1)],by.x="gene",by.y="queryID",all.x=T)
for(i in ncol(geneInf)-0:9){
        colnames(geneInf)[i] <- paste(colnames(geneInf)[i],"mOTUpresent",sep=".")
        geneInf[is.na(geneInf[,i]),i] <- "unknown"
}

#contigwise table (all contigs)
#system(paste("perl /home/users/aheintzbuschart/myScripts/calculateContigLength.pl /work/projects/ecosystem_biology/MUST/CombinedAssembly/Annotation/hmms/",FAM,"/",LIB,"/contigs.formatted.fa > /work/projects/ecosystem_biology/MUST/CombinedAssembly/Annotation/sequenceInfo/",FAM,"/",LIB,"/contigs.length.tsv",sep=""))
contigInf <- read.delim(paste("/work/projects/ecosystem_biology/MUST/CombinedAssembly/Annotation/sequenceInfo/",FAM,"/",LIB,"/contigs.length.tsv",sep=""),stringsAsFactors=F) #this contains info on ALL contigs
# DNA coverage of contigs:
cov <- read.delim(paste("/work/projects/ecosystem_biology/MUST/CombinedAssembly/Annotation/mapping/",FAM,"/",LIB,"/",dLIB,".DNAonContigs.cov.tsv",sep=""),stringsAsFactors=F,colClasses=c("character","NULL","numeric","NULL"))
contigInf <- merge(contigInf,cov,by.x=1,by.y=1,all.x=T)
colPrime <- paste("aveCov",substr(dLIB,regexpr("V",dLIB)[1],regexpr("V",dLIB)[1]+1),sep="")
colnames(contigInf)[ncol(contigInf)] <- colPrime
contigInf[is.na(contigInf[,colPrime]),colPrime] <- 0
# contig coverage by other libraries
others <- grep(dLIB,unlist(strsplit(grep(dLIB,visits$V1,value=T),split=" ")),invert=T,value=T)
for(oLIB in others){
colOth <- paste("aveCov",substr(oLIB,regexpr("V",oLIB)[1],regexpr("V",oLIB)[1]+1),sep="")
cov <- read.delim(paste("/work/projects/ecosystem_biology/MUST/CombinedAssembly/Annotation/mapping/",FAM,"/",LIB,"/",oLIB,".DNAonContigs.cov.tsv",sep=""),stringsAsFactors=F,colClasses=c("character","NULL","numeric","NULL"))
contigInf <- merge(contigInf,cov,by.x=1,by.y=1,all.x=T)
colnames(contigInf)[ncol(contigInf)] <- colOth
contigInf[is.na(contigInf[,colOth]),colOth] <- 0
}
#variants
cv <- aggregate(var$pos,list(var$contig),length)
vcv <- aggregate(var$DNAvarCov,list(var$contig),mean)
ccv <- aggregate(var$DNACov,list(var$contig),mean)
vcv <- merge(vcv,ccv,by.x=1,by.y=1)
vcv$rcv <- vcv[,2]/vcv[,3]
cv <- merge(cv,vcv[,c(1,4)],by.x=1,by.y=1)
contigInf <- merge(contigInf,cv,by.x="sequenceID",by.y=1,all.x=T)
colnames(contigInf)[ncol(contigInf)-1:0] <- c("varPerMB","varRelCov")
for(i in ncol(contigInf)-0:1){
        contigInf[is.na(contigInf[,i]),i] <- 0
}
contigInf$varPerMB <- contigInf$varPerMB*1000000/contigInf$length
pv <- aggregate(var$pos,list(var$contig),function(x) paste(unique(unlist(list(x))),sep=";",collapse=";"))
contigInf <- merge(contigInf,pv,by.x="sequenceID",by.y=1,all.x=T)
colnames(contigInf)[ncol(contigInf)] <- "varPos"
contigInf$varPos[is.na(contigInf$varPos)] <- "none"

# Kraken data:
krak <- read.delim(paste("/work/projects/ecosystem_biology/MUST/CombinedAssembly/Annotation/kraken/",FAM,"/",LIB,"/",LIB,".contigs.kraken_BacVirGenBank40GBannoUnambig.tsv",sep=""),stringsAsFactors=F)
krak <- krak[,-c(2:9)]
colnames(krak)[-1] <- paste(colnames(krak)[-1],".Kraken",sep="")
contigInf <- merge(contigInf,krak,by.x=1,by.y=1,all.x=T)
contigInf[,ncol(contigInf)-0:6][is.na(contigInf[,ncol(contigInf)-0:6])] <- "unknown" 
contigInf$annotationLevel.Kraken[is.na(contigInf$annotationLevel.Kraken)] <- "none"

#collapse gene annotations on contigs
for(anno in colnames(geneInf)[ncol(geneInf)-33:0]){
write.table(geneInf[!(geneInf[,anno] %in% c("notEssential","unknown","none","")),c("contig",anno)],paste(anno,"ByGenes.tsv",sep=""),row.names=F,quote=F,sep="\t")
system(paste("python /home/users/aheintzbuschart/myScripts/141116_MUST_gene2contig.py ",anno,"ByGenes.tsv",sep=""))
geneAnno <- read.delim(paste(anno,"ByGenes.contigs.tsv",sep=""),stringsAsFactors=F)
contigInf <- merge(contigInf,geneAnno,by.x=1,by.y=1,all.x=T)
contigInf[is.na(contigInf[,anno]),anno] <- "unknown"
}
contigInf$essentialGene[contigInf$essentialGene=="unknown"] <- "notEssential"

colnames(contigInf)[1:2] <- c("contig","length")

contigCoord <- read.delim("contigs.1000.rRNAcut.fa_5mer_clr.coords",header=F,stringsAsFactors=F,sep=",") 
contigNames <- read.delim("contigs.1000.rRNAcut.names.txt",header=F,stringsAsFactors=F)
contigNames <- sapply(gsub(">","",contigNames[,1],fixed=T),function(x)unlist(strsplit(x,split=" "))[1])
rownames(contigCoord) <- contigNames
colnames(contigCoord) <- c("x","y")

contigInf <- merge(contigInf,contigCoord,by.x=1,by.y=0,all.x=T)
contigInfo <- contigInf[!is.na(contigInf$x),]

geneInf <- merge(geneInf,contigInf,by.x="contig",by.y=1,all.x=T,suffixes=c("","contigs"))
geneInfo <- geneInf[!is.na(geneInf$x),]

rm(list=c("amp","anno","args","ccv","cv","colOth","colPrime","contigCoord","contigNames","cov","dupGenes","dupliEss","dupliID","essGenes","extrCont","i","ids","geneAnno","gvar","kegg","krak","m","metCy","mge","mphFiles","oLIB","others","pfam","protInf","pv","rnaInf","rnafile","sw","tigr","visits","var","vcv"))

save.image("WSvar.Rdata")
