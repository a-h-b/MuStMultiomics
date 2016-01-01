funcTab <- matrix(0,nrow=12,ncol=41,dimnames=list(c("geneNo","KOgenes","metaCycGenes","swissProtGenes","pfamGenes","tigrPfamGenes","annoGenes","KOrich","metaCycRich","swissProtRich","pfamRich","tigrPfamRich"),c(1:41)))
i <- 1
for(fam in paste("must_m_0",1:4,sep="")){
setwd(paste("/work/projects/ecosystem_biology/MUST/CombinedAssembly/Annotation/Bmaps/",fam,sep=""))
ids <- read.delim("ids",header=F,stringsAsFactors=F)
for(lib in ids$V2){
setwd(paste("/work/projects/ecosystem_biology/MUST/CombinedAssembly/Annotation/Bmaps/",fam,"/",lib,sep=""))
load("WSvar.Rdata")
colnames(funcTab)[i] <- lib
funcTab[1,i] <- nrow(geneInf)
funcTab[2,i] <- length(geneInf$KO[geneInf$KO!="unknown"])
funcTab[3,i] <- length(geneInf$metaCycID[geneInf$metaCycID!="unknown"])
funcTab[4,i] <- length(geneInf$swissprotEC[geneInf$swissprotEC!="unknown"])
funcTab[5,i] <- length(geneInf$pfamID[geneInf$pfamID!="unknown"])
funcTab[6,i] <- length(geneInf$tigrID[geneInf$tigrID!="unknown"])
funcTab[7,i] <- length(geneInf$KO[geneInf$KO!="unknown"|geneInf$metaCycID!="unknown"|geneInf$swissprotEC!="unknown"|geneInf$pfamID!="unknown"|geneInf$tigrID!="unknown"])
funcTab[8,i] <- length(unique(unlist(sapply(geneInf$KO[geneInf$KO!="unknown"],function(x)unlist(strsplit(x,split=";"))))))
funcTab[9,i] <- length(unique(unlist(sapply(geneInf$metaCycID[geneInf$metaCycID!="unknown"],function(x)unlist(strsplit(x,split=";"))))))
funcTab[10,i] <- length(unique(unlist(sapply(geneInf$swissprotEC[geneInf$swissprotEC!="unknown"],function(x)unlist(strsplit(x,split=";"))))))
funcTab[11,i] <- length(unique(unlist(sapply(geneInf$pfamID[geneInf$pfamID!="unknown"],function(x)unlist(strsplit(x,split=";"))))))
funcTab[12,i] <- length(unique(unlist(sapply(geneInf$tigrID[geneInf$tigrID!="unknown"],function(x)unlist(strsplit(x,split=";"))))))
i<- i+1
}
}
write.table(funcTab,"/work/projects/ecosystem_biology/MUST/CombinedAssembly/Annotation/hmms/compareAnnoStats.tsv",sep="\t",quote=F)


