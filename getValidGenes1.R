# this script is used in the MuSt project to find the names of the genes (of one COG) in the mOTU padded database which represent the mOTUs present in a sample

#inputs: as argument in the call - COG, the COG ID of the genes of interest
#       uses - the mapping files "mOTU.v1.padded.motu.linkage.map" and "mOTU.v1.map.txt" of the mOTU database and
#            - a prepared file with the names of the present mOTUs "mOTUs.present.all"  
#output: a list of genenames (written to a file named with the COG ID followed by ".present.mOTU.genenames.txt")
# 

#written by Anna Heintz-Buschart (February 2015)

args<-commandArgs(TRUE)
COG <- args[1]

phyl <- read.delim("../../mOTU.v1.padded.motu.linkage.map",stringsAsFactors=F)
map <- read.delim("../../mOTU.v1.map.txt",stringsAsFactors=F)
pres <- read.delim("mOTUs.present.all",header=F,stringsAsFactors=F)

phyl2 <- phyl[ phyl$mOTU.species.annotation %in% pres$V1,]
mp <- merge(map,phyl2,by.x="mOTU",by.y=1)
coord <- read.delim(paste("../../",COG,"/",COG,".mOTU.v1.padded.coord",sep=""),header=F,stringsAsFactors=F)
co <- gsub(" ",".",coord$V1)
coord <- read.delim(paste("../../",COG,"/",COG,".mOTU.v1.padded.coord",sep=""),header=F,stringsAsFactors=F,sep=" ")
coord$name <- co
mp2 <- merge(mp[mp$MarkerGene ==COG,],coord[,c(1,4)],by.x="SequenceID",by.y=1)

write.table(mp2$name,paste(COG,".present.mOTU.genenames.txt",sep=""),col.names=F,quote=F,sep="\t",row.names=F)



