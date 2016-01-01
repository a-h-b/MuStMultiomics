# this is the script used in the MuSt project for clustering of contigs in an automated way
# in contrast to later versions, it uses an R-workspace (WSvar.Rdata) as input which contains all the necessary information 

#outputs:
# intermittent documentation of the iterations:
#     - plots with BH-SNE maps of the clusters (on standard out)
#     - a table (reachabilityDistanceEstimates.tsv) with the used EST values
#     - a table (clusterFirstScan.tsv) with the clusters of the first iteration
#     - a table each (blob1ClusterXX.tsv) with the subclusters from the second iteration
#     - a table each (blob2ClusterXX.tsv) with the subclusters from the third iteration
#     - a table each (blob3ClusterXX.tsv) with the subclusters from the fourth iteration
#     - a table each (blob4ClusterXX.tsv) with the subclusters from the fifth iteration
#     - a table (bimodalClusterCutoffs.tsv) with the cut-off values used in separating contigs based on metagenomic coverage
# a plot with the BH-SNE map with the final cluster membership (on standard out)
# a table with the cluster membership of every contig (contigs2clusters.tsv)
# a table with stats for each cluster (clusterStats.tsv)
# a pdf file each for every population with more than 67% complete genes with a visual evaluation of the reconstructed genome
# for every cluster, a bed file each for the contigs and the genes
# R-workspaces:
#     - WSclusterAll.Rdata with all the data that was in the workspace before + the cluster information
#     - WSclusterLong.Rdata with the data on contigs >= 1000 nt that was in the workspace before + the cluster information
#     - WSclusterCompletish.Rdata with the cluster information

#written by Anna Heintz-Buschart (March 2015)

#.libPaths("/home/users/aheintzbuschart/lib/Rlibs")
library(caTools)
library(fpc)
library(FNN)
library(RColorBrewer)
library(scales)
library(diptest)
library(mixtools)

load("WSvar.Rdata")

find.cutoff <- function(data,k=2,proba=0.5) {
  model <- normalmixEM(x=data, k=k,maxrestarts=50)
  if(model!=-1){
  i <- which.min(model$mu)
  f <- function(x) {
    proba - (model$lambda[i]*dnorm(x, model$mu[i], model$sigma[i]) /(model$lambda[1]*dnorm(x, model$mu[1], model$sigma[1]) + model$lambda[2]*dnorm(x, model$mu[2], model$sigma[2])))
  }
  lower <- min(model$mu)
  upper <- max(model$mu)
  if(f(lower)*f(upper)>0){
    return(median(data))
  }else{
    return(uniroot(f=f, lower=lower, upper=upper)$root)  # Careful with division by zero if changing lower and upper
  }
}else{
  return(median(data))
}
}
evalSubPop <- function(popOIcontigs,popOIgenes,name="population",contigInf=contigInfo,geneInf=geneInfo){
evals <- list()
tableName <- paste(LIB,"_",name,"_evaluationResults.txt",sep="")
figureName <- paste(LIB,"_",name,"_evaluationResults.pdf",sep="")
evals$contigCount <- nrow(popOIcontigs)
evals$totalLength <- sum(popOIcontigs$length)
  
  pdf(figureName,width=3.5,height=4.5,pointsize=8)
  par(mar=rep(0.5,4),pty="s")
  plotlim <- 1.05*max(abs(c(min(contigInf[,ncol(contigInf)-1]),min(contigInf[,ncol(contigInf)]),max(contigInf[,ncol(contigInf)-1]),max(contigInf[,ncol(contigInf)]))))
  plot(contigInf[,ncol(contigInf)-1:0],pch=16,cex=0.3,las=1,col="black",ylim=c(-plotlim,plotlim),xlim=c(-plotlim,plotlim),ann=F,axes=F)
  points(popOIcontigs[,ncol(popOIcontigs)-1:0],pch=16,cex=0.3,las=1,col="deeppink")
  box(lwd=2)
  mtext(name,3,-1)
  par(mar=c(3,3,2.5,4),pty="m")
  d <- hist(log10(popOIcontigs$length),breaks=seq(3,6.5,length.out=50),plot=F)
  plot(d$mids, d$counts,ann=F,axes=F,xlim=c(3,6.5),ylim=c(0,1.5*max(d$counts)),col="red",type="l")
  axis(2,las=1,cex.axis=0.8,line=0,mgp=c(2,0.6,0),col.axis="red")
  par(new=T)
  d <- hist(log10(contigInf$length),breaks=seq(3,6.5,length.out=50),plot=F)
  plot(d$mids, d$counts,ann=F,axes=F,xlim=c(3,6.5),ylim=c(0,10000),type="l",lty=2)
  axis(1,at=3:6,labels=format(10^(3:6),scientific=T),cex.axis=0.8,line=0,mgp=c(2,0.6,0))
  axis(4,las=1,cex.axis=0.8,line=0,mgp=c(2,0.6,0))
  mtext("length of contigs [bp]", 1, 1.6)
  mtext("counts in population", 2, 2,col="red")
  mtext("counts in assembly", 4, 3)
  mtext("Contig length distribution",3,1,font=2)
  mtext(paste("Total length",evals$totalLength,"nt"),3,0,col="red",cex=0.9)
  box(bty="u",lty=1.5)
  
  evals$aveGC <- weighted.mean(popOIcontigs$GCperc,popOIcontigs$length)
  spec <- colorRampPalette(brewer.pal(9,"Spectral")[9:1])
  phylcol <- spec(length(unique(popOIcontigs$phylum.Kraken[popOIcontigs$phylum.Kraken!="unknown"])))[as.numeric(as.factor(popOIcontigs$phylum.Kraken[popOIcontigs$phylum.Kraken!="unknown"]))]
  par(mar=c(3,3,2,0.6),mgp=c(1.8,0.6,0))
  plot(popOIcontigs$length,popOIcontigs$GCperc*100,log="x",las=1,col="grey80",pch=16,ylim=c(5,75),bty="l",xlab="contig length",ylab="%G+C",main="Length and %G+C",cex.main=1)
  points(popOIcontigs$length[popOIcontigs$phylum.Kraken!="unknown"],popOIcontigs$GCperc[popOIcontigs$phylum.Kraken!="unknown"]*100,las=1,col=alpha(phylcol,0.8),pch=16)
  points(popOIcontigs$length[popOIcontigs$essentialGene != "notEssential"],(popOIcontigs$GCperc*100)[popOIcontigs$essentialGene != "notEssential"],pch=4,col=colors()[525])
  legend("topright",c("essential gene(s)"),title="contigs with:",bty="n",cex=0.7,y.intersp=0.8,pch=4,col=colors()[525])  
  legend("bottom",legend=c(unique(popOIcontigs$phylum.Kraken[popOIcontigs$phylum.Kraken!="unknown"]),"unknown"),text.col=alpha("black",0),col="grey80",pch=16,cex=0.7,bty="n",y.intersp=0.8,ncol=3)
  legend("bottom",legend=c(unique(popOIcontigs$phylum.Kraken[popOIcontigs$phylum.Kraken!="unknown"]),"unknown"),col=alpha(c(unique(phylcol),"grey80"),0.8),pch=16,cex=0.7,bty="n",y.intersp=0.8,ncol=3)
  colnames(popOIcontigs)[grep("aveCov",colnames(popOIcontigs))[1]] <- "aveCov"
  evals$aveCov <- weighted.mean(popOIcontigs$aveCov,popOIcontigs$length)
  par(mar=c(3,4.5,2,0.6))
	plot(popOIcontigs$length,popOIcontigs$aveCov,col="grey80",ylim=c(0.1*min(popOIcontigs$aveCov),10*max(popOIcontigs$aveCov)),log="xy",las=1,pch=16,bty="l",xlab="contig length",ylab="",main="Length and coverage",cex.main=1)
  points(popOIcontigs$length[popOIcontigs$phylum.Kraken!="unknown"],popOIcontigs$aveCov[popOIcontigs$phylum.Kraken!="unknown"],las=1,col=alpha(phylcol,0.8),pch=16)
  points(popOIcontigs$length[popOIcontigs$essentialGene != "notEssential"],popOIcontigs$aveCov[popOIcontigs$essentialGene != "notEssential"],pch=4,col=colors()[525])
  mtext("average coverage",2,3.5)
  legend("topright",c("essential gene(s)"),title="contigs with:",bty="n",cex=0.7,y.intersp=0.8,pch=4,col=colors()[525])  
  legend("bottomleft",legend=c(unique(popOIcontigs$phylum.Kraken[popOIcontigs$phylum.Kraken!="unknown"]),"unknown"),text.col=alpha("black",0),
         col="grey80",pch=16,cex=0.7,bty="n",y.intersp=0.8,ncol=3)
  legend("bottomleft",legend=c(unique(popOIcontigs$phylum.Kraken[popOIcontigs$phylum.Kraken!="unknown"]),"unknown"),
         col=alpha(c(unique(phylcol),"grey80"),0.8),pch=16,cex=0.7,bty="n",y.intersp=0.8,ncol=3) 
  evals$aveVarPerMB <- weighted.mean(popOIcontigs$varPerMB,popOIcontigs$length)
  par(mar=c(3,4.5,2,0.6))
  plot(popOIcontigs$length,popOIcontigs$varPerMB,col="grey80",ylim=c(0.1*min(popOIcontigs$varPerMB[popOIcontigs$varPerMB>0]),10*max(popOIcontigs$varPerMB[popOIcontigs$varPerMB>0])),log="xy",las=1,pch=16,bty="l",xlab="contig length",ylab="",main="Length and variant density",cex.main=1)
  points(popOIcontigs$length[popOIcontigs$phylum.Kraken!="unknown"],popOIcontigs$varPerMB[popOIcontigs$phylum.Kraken!="unknown"],las=1,col=alpha(phylcol,0.8),pch=16)
  points(popOIcontigs$length[popOIcontigs$essentialGene != "notEssential"],popOIcontigs$varPerMB[popOIcontigs$essentialGene != "notEssential"],pch=4,col=colors()[525])
  mtext("variant density",2,3.5)
  legend("topright",c("essential gene(s)"),title="contigs with:",bty="n",cex=0.7,y.intersp=0.8,pch=4,col=colors()[525])
  legend("bottomleft",legend=c(unique(popOIcontigs$phylum.Kraken[popOIcontigs$phylum.Kraken!="unknown"]),"unknown"),text.col=alpha("black",0),col="grey80",pch=16,cex=0.7,bty="n",y.intersp=0.8,ncol=3)
  legend("bottomleft",legend=c(unique(popOIcontigs$phylum.Kraken[popOIcontigs$phylum.Kraken!="unknown"]),"unknown"),col=alpha(c(unique(phylcol),"grey80"),0.8),pch=16,cex=0.7,bty="n",y.intersp=0.8,ncol=3) 
  
  par(mar=c(3,3,2.5,4))
  d <- hist(log10(popOIcontigs$aveCov[popOIcontigs$aveCov>=0.1]),breaks=seq(-1,3.8,length.out=50),plot=F)
  plot(d$mids, d$counts,ann=F,axes=F,xlim=c(-1,3.8),ylim=c(0,1.5*max(d$counts)),col="red",type="l")
  axis(2,las=1,cex.axis=0.8,line=0,mgp=c(2,0.6,0),col.axis="red")
  par(new=T)
  d <- hist(log10(contigInf[contigInf[,4]>=0.1&contigInf[,4]<=1000,4]),breaks=seq(-1,3.8,length.out=50),plot=F)
  plot(d$mids, d$counts,ann=F,axes=F,xlim=c(-1,3.8),ylim=c(0,10000),type="l",lty=2)
  axis(1,at=-1:3,labels=format(10^(-1:3),scientific=T),cex.axis=0.8,line=0,mgp=c(2,0.6,0))
  axis(4,las=1,cex.axis=0.8,line=0,mgp=c(2,0.6,0))
  mtext("coverage of contigs", 1, 2)
  mtext("counts in population", 2, 2,col="red")
  mtext("counts in assembly", 4, 3)
  mtext("Contig coverage distribution",3,1,font=2)
  mtext(paste(evals$contigCount,"contigs; ",length(which(popOIcontigs$aveCov<0.1)),"contigs without coverage"),3,0,col="red",cex=0.9)
  box(bty="u",lty=1.5)
  
  par(mar=c(3,3,2.5,4))
  d <- hist(log10(popOIcontigs$varPerMB[popOIcontigs$varPerMB>=1]),breaks=seq(0,5.5,length.out=50),plot=F)
  plot(d$mids, d$counts,ann=F,axes=F,xlim=c(0,5.5),ylim=c(0,1.5*max(d$counts)),col="red",type="l")
  axis(2,las=1,cex.axis=0.8,line=0,mgp=c(2,0.6,0),col.axis="red")
  par(new=T)
  d <- hist(log10(contigInf$varPerMB[contigInf$varPerMB>=1&contigInf$varPerMB<10^5.5]),breaks=seq(0,5.5,length.out=50),plot=F)
  plot(d$mids, d$counts,ann=F,axes=F,xlim=c(0,5.5),ylim=c(0,10000),type="l",lty=2)
  axis(1,at=0:5,labels=format(10^(0:5),scientific=T),cex.axis=0.8,line=0,mgp=c(2,0.6,0))
  axis(4,las=1,cex.axis=0.8,line=0,mgp=c(2,0.6,0))
  mtext("variant density of contigs", 1, 2)
  mtext("counts in population", 2, 2,col="red")
  mtext("counts in assembly", 4, 3)
  mtext("Contig variant density distribution",3,1,font=2)
  mtext(paste(evals$contigCount,"contigs; ",length(which(popOIcontigs$varPerMB<1)),"contigs without variants"),3,0,col="red",cex=0.9)
  box(bty="u",lty=1.5)
  
  evals$geneCount <- nrow(popOIgenes)
  evals$completeGeneCount <- length(which(popOIgenes$completeness=="complete"))  
  
  par(mar=c(3,3,2.5,4))
  d <- hist(log10(popOIgenes$length),breaks=seq(1,5,length.out=50),plot=F)
  helpd <- d$counts
  plot(d$mids, helpd,ann=F,axes=F,xlim=c(1,5),ylim=c(0,1.5*max(helpd)),col="red",type="l")
  axis(2,las=1,cex.axis=0.8,line=0,mgp=c(2,0.6,0),col.axis="red")
  par(new=T)
  d <- hist(log10(popOIgenes$length[ popOIgenes$completeness=="complete"]),breaks=seq(1,5,length.out=50),plot=F)
  plot(d$mids, d$counts,ann=F,axes=F,xlim=c(1,5),ylim=c(0,1.5*max(helpd)),col="red",type="l",lty=2)
  axis(2,las=1,cex.axis=0.8,line=0,mgp=c(2,0.6,0),col.axis="red")
  par(new=T)
  d <- hist(log10(geneInf$length),breaks=seq(1,5,length.out=50),plot=F)
  plot(d$mids, d$counts,ann=F,axes=F,xlim=c(1,5),ylim=c(0,35000),type="l",lty=1)
  par(new=T)
  d <- hist(log10(geneInf$length[geneInf$completeness=="complete"]),breaks=seq(1,5,length.out=50),plot=F)
  plot(d$mids, d$counts,ann=F,axes=F,xlim=c(1,5),ylim=c(0,35000),type="l",lty=2)
  axis(1,at=1:5,labels=format(10^(1:5),scientific=T),cex.axis=0.8,line=0,mgp=c(2,0.5,0))
  axis(4,las=1,cex.axis=0.8,line=0,mgp=c(2,0.5,0))
  mtext("length of genes [bp]", 1, 2)
  mtext("counts in population", 2, 2,col="red")
  mtext("counts in assembly", 4, 3)
  mtext("Gene length distribution",3,1,font=2)
  box(bty="u",lty=1.5)
  legend("topleft",c("complete","incomplete"),lty=c(1,2),cex=0.8,bty="n",col="red")
  legend("topright",c("complete","incomplete"),lty=c(1,2),cex=0.8,bty="n")
  
  par(mar=c(5,5,5,5))
  pie(c(evals$completeGeneCount,evals$geneCount-evals$completeGeneCount),labels=paste(c("complete","incomplete")," (",c(evals$completeGeneCount,evals$geneCount-evals$completeGeneCount),")",sep=""),
      col=c("grey20","grey80"),clockwise=T,main="Completeness of genes")
  
  evals$expressedGeneCount <- length(which(popOIgenes$aveCovRNAfw>0))
  evals$antisenseGeneCount <- length(which(popOIgenes$aveCovRNArc>0)) 
  par(mar=c(3,5,2,0.6),mgp=c(3.8,0.6,0))
  plot(ifelse(popOIgenes$aveCovRNAfw>1,log10(popOIgenes$aveCovRNAfw),0),type="h",xlab="",
       ylab="average coverage RNA",main="Gene coverage with RNA reads",cex.main=1.2,axes=F,ylim=c(-4.5,4.5))
  points((1:length(popOIgenes$aveCovRNArc)),
         ifelse(popOIgenes$aveCovRNArc>1,-log10(popOIgenes$aveCovRNArc),0),type="h")
  mtext("gene index",1,1.7)
  mtext("antisense",2,2.8,adj=0.2)
  mtext("sense",2,2.8,adj=0.8)
  box(bty="l",lwd=1.5)
  axis(1)
  axis(2,at=c(-4,-2,0,2,4),las=1,labels=format(10^c(4,2,0,2,4),scientific=T))
  
  evals$rRNACount <- length(which(popOIgenes$kind!="protein"))
  evals$KOGeneCount <- length(which(popOIgenes$KO!="unknown"))
  evals$metaCycCount <- length(which(popOIgenes$metaCycID!="unknown"))
  evals$swissprotCount <- length(which(popOIgenes$swissprotEC!="unknown"))
  evals$pfamCount <- length(which(popOIgenes$pfamID!="unknown"))
  evals$tigrCount <- length(which(popOIgenes$tigrID!="unknown"))
  evals$annoCount <- length(which(popOIgenes$tigrID!="unknown"|popOIgenes$pfamID!="unknown"|popOIgenes$swissprotEC!="unknown"|popOIgenes$metaCycID!="unknown"|popOIgenes$KO!="unknown"))
  par(mar=c(5,5,5,5))
  pie(c(evals$KOGeneCount,evals$geneCount-evals$KOGeneCount),labels=paste(c("annotated","unknown")," (",c(evals$KOGeneCount,evals$geneCount-evals$KOGeneCount),")",sep=""),
      col=c("grey20","grey80"),clockwise=T,main="KO-annotation of genes")
  pie(c(evals$metaCycCount,evals$geneCount-evals$metaCycCount),labels=paste(c("annotated","unknown")," (",c(evals$metaCycCount,evals$geneCount-evals$metaCycCount),")",sep=""),
      col=c("grey20","grey80"),clockwise=T,main="MetaCyc-annotation of genes")
  pie(c(evals$swissprotCount,evals$geneCount-evals$swissprotCount),labels=paste(c("annotated","unknown")," (",c(evals$swissprotCount,evals$geneCount-evals$swissprotCount),")",sep=""),
      col=c("grey20","grey80"),clockwise=T,main="Swiss-Prot-annotation of genes")
  pie(c(evals$pfamCount,evals$geneCount-evals$pfamCount),labels=paste(c("annotated","unknown")," (",c(evals$pfamCount,evals$geneCount-evals$pfamCount),")",sep=""),
      col=c("grey20","grey80"),clockwise=T,main="Pfam-annotation of genes")
  pie(c(evals$tigrCount,evals$geneCount-evals$tigrCount),labels=paste(c("annotated","unknown")," (",c(evals$tigrCount,evals$geneCount-evals$tigrCount),")",sep=""),
      col=c("grey20","grey80"),clockwise=T,main="tigrPfam-annotation of genes")
  pie(c(evals$annoCount,evals$geneCount-evals$annoCount),labels=paste(c("annotated","unknown")," (",c(evals$annoCount,evals$geneCount-evals$annoCount),")",sep=""),
      col=c("grey20","grey80"),clockwise=T,main="any annotation of genes")
  
  write.table(t(c("Contig Count","Mean Contig Length","Ave Contig %GC", "Ave Contig Cov")),tableName,row.names=F,col.names=F,quote=F,sep="\t")
  write.table(t(c(evals$contigCount,evals$totalLength,evals$aveGC,evals$aveCov)),tableName,row.names=F,col.names=F,quote=F,append=T,sep="\t")
  write.table(t(c("Gene Count","rRNA Gene Count", "Complete Gene Count","Expressed Gene Count","Antisense-expressed Gene Count")),tableName,row.names=F,col.names=F,quote=F,sep="\t")
  write.table(t(c(evals$geneCount,evals$rRNACount,evals$completeGeneCount,evals$expressedGeneCount,evals$antisenseGeneCount)),tableName,row.names=F,col.names=F,quote=F,append=T,sep="\t")
  write.table(t(c("annotated Gene Count","KO-annotated Gene Count","MetaCyc-annotated Gene Count","Swiss-Prot-annotated Gene Count","Pfam-annotated Gene Count","Tigr-Pfam-annotated Gene Count")),tableName,row.names=F,col.names=F,quote=F,sep="\t")
  write.table(t(c(evals$anooCount,evals$KOGeneCount,evals$metaCycCount,evals$swissprotCount,evals$pfamCount,evals$tigrCount)),tableName,row.names=F,col.names=F,quote=F,append=T,sep="\t")
  
  evals$essentialGeneTable <- table(table(popOIgenes$essentialGene[popOIgenes$essentialGene!="notEssential"]))
  pie(c(evals$essentialGeneTable,109-sum(evals$essentialGeneTable)),labels=c(names(evals$essentialGeneTable),"missing"),
      col=brewer.pal(length(evals$essentialGeneTable)+1,"Set1"),clockwise=T,main="Presence of essential genes")
  
  write.table("Essential Genes for Completeness Analysis",tableName,row.names=F,col.names=F,quote=F,append=T,sep="\t")
  write.table(rbind(names(evals$essentialGeneTable),t(evals$essentialGeneTable)),tableName,row.names=F,col.names=F,quote=F,append=T,sep="\t")
  
  
  evals$motuGeneTable <- table(table(popOIgenes$MarkerGene.mOTUpresent[popOIgenes$MarkerGene.mOTUpresent!="unknown"]))
  write.table("Marker Genes linked to mOTUs",tableName,row.names=F,col.names=F,quote=F,append=T,sep="\t")
  write.table(rbind(names(evals$motuGeneTable),t(evals$motuGeneTable)),tableName,row.names=F,col.names=F,quote=F,append=T,sep="\t")
  
  moKing <- popOIgenes$Superkingdom.mOTUpresent[ popOIgenes$MarkerGene.mOTUpresent!="unknown"]
  if(length(moKing)>0){
    evals$kingdomMotuTab <- sort(table(moKing),decreasing=T)
    
    write.table("Superkingdom-level Annotation of mOTU Genes",tableName,row.names=F,col.names=F,quote=F,append=T,sep="\t")
    write.table(rbind(names(evals$kingdomMotuTab),t(evals$kingdomMotuTab)),tableName,row.names=F,col.names=F,quote=F,append=T,sep="\t")
    
    plab <- paste(names(evals$kingdomMotuTab)," (",evals$kingdomMotuTab,")",sep="")
    par(mar=c(5,5,5,5))
    pie(evals$kingdomMotuTab,labels=plab,col=brewer.pal(9,"Set1"),clockwise=T)
    mtext("Superkingdom-level annotation (mOTUs)",3,1,font=2,cex=1.2)

    moPhy <- popOIgenes$Phylum.mOTUpresent[ popOIgenes$MarkerGene.mOTUpresent!="unknown"]
    evals$phylumMotuTab <- sort(table(moPhy),decreasing=T)
    
    write.table("Phylum-level Annotation of mOTU Genes",tableName,row.names=F,col.names=F,quote=F,append=T,sep="\t")
    write.table(rbind(names(evals$phylumMotuTab),t(evals$phylumMotuTab)),tableName,row.names=F,col.names=F,quote=F,append=T,sep="\t")
    
    trCol <- data.frame("phylum"=unique(popOIcontigs$phylum.Kraken[popOIcontigs$phylum.Kraken!="unknown"]),
                        "col"=unique(phylcol),stringsAsFactors=F)
    trCol <- trCol[trCol$phylum %in% unique(moPhy),]
    trCol <- rbind(trCol,c("unknown","grey80"))
    colnames(trCol) <- c("phylum","col")
    addPhyl <- setdiff(unique(moPhy),
                       unique(popOIcontigs$phylum.Kraken[popOIcontigs$phylum.Kraken!="unknown"]))
    if(length(addPhyl)>0){
      addSpec <- colorRampPalette(brewer.pal(9,"RdPu")[7:1])
      addCol <- addSpec(length(addPhyl))
      trCol <- rbind(trCol,data.frame("phylum"=addPhyl,"col"=addCol,stringsAsFactors=F))
    }
    plab <- paste(names(evals$phylumMotuTab)," (",evals$phylumMotuTab,")",sep="")
    par(mar=c(5,5,5,5))
    pie(evals$phylumMotuTab,labels=plab,
        col=unlist(sapply(names(evals$phylumMotuTab),function(x) trCol$col[trCol$phylum==x])),clockwise=T)
    mtext("Phylum level annotation (mOTUs)",3,1,font=2)
    
    moClass <- popOIgenes$Class.mOTUpresent[ popOIgenes$MarkerGene.mOTUpresent!="unknown"]
    evals$classMotuTab <- sort(table(moClass),decreasing=T)
    
    write.table("Class-level Annotation of mOTU Genes",tableName,row.names=F,col.names=F,quote=F,append=T,sep="\t")
    write.table(rbind(names(evals$classMotuTab),t(evals$classMotuTab)),tableName,row.names=F,col.names=F,quote=F,append=T,sep="\t")
    
    plab <- paste(names(evals$classMotuTab)," (",evals$classMotuTab,")",sep="")
    par(mar=c(6,6,6,6))
    pie(evals$classMotuTab,labels=plab,col=brewer.pal(9,"Set1"),clockwise=T)
    mtext("Class level annotation (mOTUs)",3,1,font=2)

    moOrd <- popOIgenes$Order.mOTUpresent[ popOIgenes$MarkerGene.mOTUpresent!="unknown"]
    evals$orderMotuTab <- sort(table(moOrd),decreasing=T)
    
    write.table("Order-level Annotation of mOTU Genes",tableName,row.names=F,col.names=F,quote=F,append=T,sep="\t")
    write.table(rbind(names(evals$orderMotuTab),t(evals$orderMotuTab)),tableName,row.names=F,col.names=F,quote=F,append=T,sep="\t")
    
    plab <- paste(names(evals$orderMotuTab)," (",evals$orderMotuTab,")",sep="")
    par(mar=c(5,5,5,5))
    pie(evals$orderMotuTab,labels=plab,col=brewer.pal(9,"Set1"),clockwise=T)
    mtext("Order level annotation (mOTUs)",3,1,font=2)
    
    moFam <- popOIgenes$Family.mOTUpresent[ popOIgenes$MarkerGene.mOTUpresent!="unknown"]
    evals$familyMotuTab <- sort(table(moFam),decreasing=T)
    
    write.table("Family-level Annotation of mOTU Genes",tableName,row.names=F,col.names=F,quote=F,append=T,sep="\t")
    write.table(rbind(names(evals$familyMotuTab),t(evals$familyMotuTab)),tableName,row.names=F,col.names=F,quote=F,append=T,sep="\t")
    
    plab <- paste(names(evals$familyMotuTab)," (",evals$familyMotuTab,")",sep="")
    par(mar=c(5,5,5,5))
    pie(evals$familyMotuTab,labels=plab,col=brewer.pal(9,"Set1"),clockwise=T)
    mtext("Family level annotation",3,1,font=2)
    
    moGenus <- popOIgenes$Genus.mOTUpresent[ popOIgenes$MarkerGene.mOTUpresent!="unknown"]
    evals$genusMotuTab <- sort(table(moGenus),decreasing=T)
    
    write.table("Genus-level Annotation of mOTU Genes",tableName,row.names=F,col.names=F,quote=F,append=T,sep="\t")
    write.table(rbind(names(evals$genusMotuTab),t(evals$genusMotuTab)),tableName,row.names=F,col.names=F,quote=F,append=T,sep="\t")
    
    plab <- paste(names(evals$genusMotuTab)," (",evals$genusMotuTab,")",sep="")
    par(mar=c(5,5,5,5))
    pie(evals$genusMotuTab,labels=plab,col=brewer.pal(9,"Set1"),clockwise=T)
    mtext("Genus level annotation (mOTUs)",3,1,font=2)
    
    moSpec <- popOIgenes$SpeciesCluster.mOTUpresent[ popOIgenes$MarkerGene.mOTUpresent!="unknown"]
    evals$speciesMotuTab <- sort(table(moSpec),decreasing=T)
    
    write.table("Species-level Annotation of mOTU Genes",tableName,row.names=F,col.names=F,quote=F,append=T,sep="\t")
    write.table(rbind(names(evals$speciesMotuTab),t(evals$speciesMotuTab)),tableName,row.names=F,col.names=F,quote=F,append=T,sep="\t")
    
    plab <- paste(names(evals$speciesMotuTab)," (",evals$speciesMotuTab,")",sep="")
    par(mar=c(5,5,5,5))
    pie(evals$speciesMotuTab,labels=plab,col=brewer.pal(9,"Set1"),clockwise=T)
    mtext("Species level annotation (mOTUs)",3,1,font=2)
                       
    motu <- popOIgenes$mOTU.species.annotation.mOTUpresent[ popOIgenes$MarkerGene.mOTUpresent!="unknown"]
    evals$motuTab <- sort(table(motu),decreasing=T)
                                          
    write.table("mOTU-level Annotation of mOTU Genes",tableName,row.names=F,col.names=F,quote=F,append=T,sep="\t")
    write.table(rbind(names(evals$motuTab),t(evals$motuTab)),tableName,row.names=F,col.names=F,quote=F,append=T,sep="\t")
    
    plab <- paste(names(evals$motuTab)," (",evals$motuTab,")",sep="")
    par(mar=c(5,5,5,5))
    pie(evals$motuTab,labels=plab,col=brewer.pal(9,"Set1"),clockwise=T)
    mtext("mOTU level annotation",3,1,font=2)
  } else {
    write.table("No Marker Genes for mOTUs in Selection",tableName,row.names=F,col.names=F,quote=F,append=T,sep="\t")
  }

  
  evals$markerGeneTable <- table(table(popOIgenes$marker.Amphora[popOIgenes$marker.Amphora!="unknown"]))
  write.table("Marker Genes for Phylogenetic Analysis",tableName,row.names=F,col.names=F,quote=F,append=T,sep="\t")
  write.table(rbind(names(evals$markerGeneTable),t(evals$markerGeneTable)),tableName,row.names=F,col.names=F,quote=F,append=T,sep="\t")
  
  amphKing <- t(sapply(popOIgenes$kingdom.Amphora[ popOIgenes$marker.Amphora!="unknown"],
                       function(x) unlist(strsplit(gsub("unknown","unknown(0",gsub("\\)","",x)),split="\\("))))
  if(ncol(amphKing)>0){
    evals$kingdomAmphoraTab <- sort(table(amphKing[,1]),decreasing=T)
    
    write.table("Kingdom-level Annotation of Marker Genes",tableName,row.names=F,col.names=F,quote=F,append=T,sep="\t")
    write.table(rbind(names(evals$kingdomAmphoraTab),t(evals$kingdomAmphoraTab)),tableName,row.names=F,col.names=F,quote=F,append=T,sep="\t")
    
    plab <- paste(names(evals$kingdomAmphoraTab)," (",evals$kingdomAmphoraTab,")",sep="")
    par(mar=c(5,5,5,5))
    pie(evals$kingdomAmphoraTab,labels=plab,col=brewer.pal(9,"Set1"),clockwise=T)
    mtext("Kingdom level annotation (amphora2)",3,1,font=2,cex=1.2)
    par(mar=c(10,4,2.5,1),mgp=c(3,0.6,0))
    boxplot(as.numeric(amphKing[,2])~amphKing[,1],las=2,names=levels(as.factor(amphKing[,1])),ylab="confidence",varwidth=T,
            main="confidence in Amphora2 annotations",border=brewer.pal(9,"Set1")[order(table(amphKing[,1]),decreasing=T)])
    if(length(unique(amphKing[,1]))==1) mtext(unique(amphKing[,1]),1,0.6,cex=0.9)
    mtext("Kingdom",1,8.8,font=2)
    amphPhy <- t(sapply(popOIgenes$phylum.Amphora[ popOIgenes$marker.Amphora!="unknown"],
                        function(x) unlist(strsplit(gsub("unknown","unknown(0",gsub("\\)","",x)),split="\\("))))
    evals$phylumAmphoraTab <- sort(table(amphPhy[,1]),decreasing=T)
    
    write.table("Phylum-level Annotation of Marker Genes",tableName,row.names=F,col.names=F,quote=F,append=T,sep="\t")
    write.table(rbind(names(evals$phylumAmphoraTab),t(evals$phylumAmphoraTab)),tableName,row.names=F,col.names=F,quote=F,append=T,sep="\t")
    
    trCol <- data.frame("phylum"=unique(popOIcontigs$phylum.Kraken[popOIcontigs$phylum.Kraken!="unknown"]),
                        "col"=unique(phylcol),stringsAsFactors=F)
    trCol <- trCol[trCol$phylum %in% unique(amphPhy[,1]),]
    trCol <- rbind(trCol,c("unknown","grey80"))
    addPhyl <- setdiff(unique(amphPhy[,1]),
                       unique(popOIcontigs$phylum.Kraken[popOIcontigs$phylum.Kraken!="unknown"]))
    if(length(addPhyl)>0){
      addSpec <- colorRampPalette(brewer.pal(9,"RdPu")[7:1])
      addCol <- addSpec(length(addPhyl))
      trCol <- rbind(trCol,data.frame("phylum"=addPhyl,"col"=addCol,stringsAsFactors=F))
    }
    plab <- paste(names(evals$phylumAmphoraTab)," (",evals$phylumAmphoraTab,")",sep="")
    par(mar=c(5,5,5,5))
    pie(evals$phylumAmphoraTab,labels=plab,
        col=unlist(sapply(names(evals$phylumAmphoraTab),function(x) trCol$col[trCol$phylum==x])),clockwise=T)
    mtext("Phylum level annotation (amphora2)",3,1,font=2)
    par(mar=c(10,4,2.5,1),mgp=c(3,0.6,0))
    boxplot(as.numeric(amphPhy[,2])~amphPhy[,1],las=2,names=levels(as.factor(amphPhy[,1])),ylab="confidence",varwidth=T,
            main="confidence in Amphora2 annotations",
            border=trCol[order(trCol$phylum),2])
    if(length(unique(amphPhy[,1]))==1) mtext(unique(amphPhy[,1]),1,0.6,cex=0.9)
    mtext("Phylum",1,8.8,font=2)
    
    amphClass <- t(sapply(popOIgenes$class.Amphora[ popOIgenes$marker.Amphora!="unknown"],
                          function(x) unlist(strsplit(gsub("unknown","unknown(0",gsub("\\)","",x)),split="\\("))))
    evals$classAmphoraTab <- sort(table(amphClass[,1]),decreasing=T)
    
    write.table("Class-level Annotation of Marker Genes",tableName,row.names=F,col.names=F,quote=F,append=T,sep="\t")
    write.table(rbind(names(evals$classAmphoraTab),t(evals$classAmphoraTab)),tableName,row.names=F,col.names=F,quote=F,append=T,sep="\t")
    
    plab <- paste(names(evals$classAmphoraTab)," (",evals$classAmphoraTab,")",sep="")
    par(mar=c(6,6,6,6))
    pie(evals$classAmphoraTab,labels=plab,col=brewer.pal(9,"Set1"),clockwise=T)
    mtext("Class level annotation (amphora2)",3,1,font=2)
    par(mar=c(10,4,2.5,1),mgp=c(3,0.6,0))
    boxplot(as.numeric(amphClass[,2])~amphClass[,1],las=2,names=levels(as.factor(amphClass[,1])),ylab="confidence",varwidth=T,
            main="confidence in Amphora2 annotations",border=brewer.pal(9,"Set1")[order(table(amphClass[,1]),decreasing=T)])
    if(length(unique(amphClass[,1]))==1) mtext(unique(amphClass[,1]),1,0.6,cex=0.9)
    mtext("Class",1,8.8,font=2)
    
    amphOrd <- t(sapply(popOIgenes$order.Amphora[ popOIgenes$marker.Amphora!="unknown"],
                        function(x) unlist(strsplit(gsub("unknown","unknown(0",gsub("\\)","",x)),split="\\("))))
    evals$orderAmphoraTab <- sort(table(amphOrd[,1]),decreasing=T)
    
    write.table("Order-level Annotation of Marker Genes",tableName,row.names=F,col.names=F,quote=F,append=T,sep="\t")
    write.table(rbind(names(evals$orderAmphoraTab),t(evals$orderAmphoraTab)),tableName,row.names=F,col.names=F,quote=F,append=T,sep="\t")
    
    plab <- paste(names(evals$orderAmphoraTab)," (",evals$orderAmphoraTab,")",sep="")
    par(mar=c(5,5,5,5))
    pie(evals$orderAmphoraTab,labels=plab,col=brewer.pal(9,"Set1"),clockwise=T)
    mtext("Order level annotation (amphora2)",3,1,font=2)
    par(mar=c(10,4,2.5,1),mgp=c(3,0.6,0))  
    boxplot(as.numeric(amphOrd[,2])~amphOrd[,1],las=2,names=levels(as.factor(amphOrd[,1])),ylab="confidence",varwidth=T,
            main="confidence in Amphora2 annotations",border=brewer.pal(9,"Set1")[order(table(amphOrd[,1]),decreasing=T)])
    if(length(unique(amphOrd[,1]))==1) mtext(unique(amphOrd[,1]),1,0.6,cex=0.9)
    mtext("Order",1,8.8,font=2)
    
    amphFam <- t(sapply(popOIgenes$family.Amphora[ popOIgenes$marker.Amphora!="unknown"],
                        function(x) unlist(strsplit(gsub("unknown","unknown(0",gsub("\\)","",x)),split="\\("))))
    evals$familyAmphoraTab <- sort(table(amphFam[,1]),decreasing=T)
    
    write.table("Family-level Annotation of Marker Genes",tableName,row.names=F,col.names=F,quote=F,append=T,sep="\t")
    write.table(rbind(names(evals$familyAmphoraTab),t(evals$familyAmphoraTab)),tableName,row.names=F,col.names=F,quote=F,append=T,sep="\t")
    
    plab <- paste(names(evals$familyAmphoraTab)," (",evals$familyAmphoraTab,")",sep="")
    par(mar=c(5,5,5,5))
    pie(evals$familyAmphoraTab,labels=plab,col=brewer.pal(9,"Set1"),clockwise=T)
    mtext("Family level annotation (amphora2)",3,1,font=2)
    par(mar=c(10,4,2.5,1),mgp=c(3,0.6,0))
    boxplot(as.numeric(amphFam[,2])~amphFam[,1],las=2,names=levels(as.factor(amphFam[,1])),ylab="confidence",varwidth=T,
            main="confidence in Amphora2 annotations",
            border=brewer.pal(9,"Set1")[as.numeric(unique(as.factor(amphFam[,1])))[order(names(evals$familyAmphoraTab))]])
    if(length(unique(amphFam[,1]))==1) mtext(unique(amphFam[,1]),1,0.6,cex=0.9)
    mtext("Family",1,8.8,font=2)
    
    amphGenus <- t(sapply(popOIgenes$genus.Amphora[ popOIgenes$marker.Amphora!="unknown"],
                          function(x) unlist(strsplit(gsub("unknown","unknown(0",gsub("\\)","",x)),split="\\("))))
    evals$genusAmphoraTab <- sort(table(amphGenus[,1]),decreasing=T)
    
    write.table("Genus-level Annotation of Marker Genes",tableName,row.names=F,col.names=F,quote=F,append=T,sep="\t")
    write.table(rbind(names(evals$genusAmphoraTab),t(evals$genusAmphoraTab)),tableName,row.names=F,col.names=F,quote=F,append=T,sep="\t")
    
    plab <- paste(names(evals$genusAmphoraTab)," (",evals$genusAmphoraTab,")",sep="")
    par(mar=c(5,5,5,5))
    pie(evals$genusAmphoraTab,labels=plab,col=brewer.pal(9,"Set1"),clockwise=T)
    mtext("Genus level annotation (amphora2)",3,1,font=2)
    par(mar=c(10,4,2.5,1),mgp=c(3,0.6,0))
    boxplot(as.numeric(amphGenus[,2])~amphGenus[,1],las=2,names=levels(as.factor(amphGenus[,1])),ylab="confidence",varwidth=T,
            main="confidence in Amphora2 annotations",border=brewer.pal(9,"Set1")[order(table(amphGenus[,1]),decreasing=T)])
    if(length(unique(amphGenus[,1]))==1) mtext(unique(amphGenus[,1]),1,0.6,cex=0.9)
    mtext("Genus",1,8.8,font=2)
    
    amphSpec <- t(sapply(popOIgenes$species.Amphora[ popOIgenes$marker.Amphora!="unknown"],
                         function(x) unlist(strsplit(gsub("unknown","unknown(0",gsub("\\)","",x)),split="\\("))))
    evals$speciesAmphoraTab <- sort(table(amphSpec[,1]),decreasing=T)
    
    write.table("Species-level Annotation of Marker Genes",tableName,row.names=F,col.names=F,quote=F,append=T,sep="\t")
    write.table(rbind(names(evals$speciesAmphoraTab),t(evals$speciesAmphoraTab)),tableName,row.names=F,col.names=F,quote=F,append=T,sep="\t")
    
    plab <- paste(names(evals$speciesAmphoraTab)," (",evals$speciesAmphoraTab,")",sep="")
    par(mar=c(5,5,5,5))
    pie(evals$speciesAmphoraTab,labels=plab,col=brewer.pal(9,"Set1"),clockwise=T)
    mtext("Species level annotation (amphora2)",3,1,font=2)
    par(mar=c(10,4,2.5,1),mgp=c(3,0.6,0))
    boxplot(as.numeric(amphSpec[,2])~amphSpec[,1],las=2,names=gsub(" ","\n",levels(as.factor(amphSpec[,1]))),ylab="confidence",varwidth=T,
            main="confidence in Amphora2 annotations",border=brewer.pal(9,"Set1")[order(table(amphSpec[,1]),decreasing=T)])
    if(length(unique(amphSpec[,1]))==1) mtext(unique(amphSpec[,1]),1,0.6,cex=0.9)
    mtext("Species",1,8.8,font=2)
  } else {
    write.table("No Amphora Marker Genes in Selection",tableName,row.names=F,col.names=F,quote=F,append=T,sep="\t")
  }  
  evals$kingdomTab <- sort(table(popOIcontigs$kingdom.Kraken),decreasing=T)
  
  write.table("Kingdom-level Annotation of Contigs",tableName,row.names=F,col.names=F,quote=F,append=T,sep="\t")
  write.table(rbind(names(evals$kingdomTab),t(evals$kingdomTab)),tableName,row.names=F,col.names=F,quote=F,append=T,sep="\t")
  
  kingdomTabNt <- aggregate(popOIcontigs$length,list(popOIcontigs$kingdom.Kraken),sum)
  rownames(kingdomTabNt) <- kingdomTabNt$Group.1
  evals$kingdomTabNt <- kingdomTabNt$x
  names(evals$kingdomTabNt) <- kingdomTabNt$Group.1
  evals$kingdomTabNt <- sort(evals$kingdomTabNt,decreasing=T)
  par(mfcol=c(2,1),mar=c(0.2,1,2.4,1))
  plab <- names(evals$kingdomTab)
  pie(evals$kingdomTab,labels=plab,col=brewer.pal(9,"Set1"),clockwise=T)
  mtext("Kingdom level annotation (Kraken)",3,1,font=2)
  mtext("contig number",3,0.2,cex=0.9)
  par(mar=c(1.1,1,1.5,1))
  plab <- sort(unique(popOIcontigs$kingdom.Kraken))[order(kingdomTabNt$x,decreasing=T)]
  plab <- evals$kingdomTabNt
  pie(evals$kingdomTabNt,labels=plab,col=brewer.pal(9,"Set1"),clockwise=T)
  mtext("contig length",3,0.2,cex=0.9)
  
  write.table("Summary of Lengths of Kingdom-level Annotation of Contigs",tableName,row.names=F,col.names=F,quote=F,append=T,sep="\t")
  write.table(rbind(names(evals$kingdomTabNt),t(evals$kingdomTabNt)),tableName,row.names=F,col.names=F,quote=F,append=T,sep="\t")
  
  
  evals$phylumTab <- sort(table(popOIcontigs$phylum.Kraken),decreasing=T)
  
  write.table("Phylum-level Annotation of Contigs",tableName,row.names=F,col.names=F,quote=F,append=T,sep="\t")
  write.table(rbind(names(evals$phylumTab),t(evals$phylumTab)),tableName,row.names=F,col.names=F,quote=F,append=T,sep="\t")
  
  phylumTabNt <- aggregate(popOIcontigs$length,list(popOIcontigs$phylum.Kraken),sum)
  rownames(phylumTabNt) <- phylumTabNt$Group.1
  evals$phylumTabNt <- phylumTabNt$x
  names(evals$phylumTabNt) <- phylumTabNt$Group.1
  evals$phylumTabNt <- sort(evals$phylumTabNt,decreasing=T)
  par(mar=c(0.2,1,2.4,1))
  plab <- sort(c(unique(popOIcontigs$phylum.Kraken[popOIcontigs$phylum.Kraken!="unknown"]),"unknown"))[order(table(popOIcontigs$phylum.Kraken),decreasing=T)]
  pie(evals$phylumTab,labels=plab,clockwise=T,
      col=c(unique(phylcol),"grey80")[order(c(unique(popOIcontigs$phylum.Kraken[popOIcontigs$phylum.Kraken!="unknown"]),"unknown"))][order(table(popOIcontigs$phylum.Kraken),decreasing=T)])
  mtext("Phylum level annotation (Kraken)",3,1.1,font=2,cex=1.2)
  mtext("contig number",3,0.2,cex=0.9,font=2)
  par(mar=c(1.1,1,1.5,1))
  plab <- sort(c(unique(popOIcontigs$phylum.Kraken[popOIcontigs$phylum.Kraken!="unknown"]),"unknown"))[order(phylumTabNt$x,decreasing=T)]
  pie(evals$phylumTabNt,labels=plab,clockwise=T,
      col=c(unique(phylcol),"grey80")[order(c(unique(popOIcontigs$phylum.Kraken[popOIcontigs$phylum.Kraken!="unknown"]),"unknown"))][order(phylumTabNt$x,decreasing=T)])
  mtext("contig length",3,0.2,cex=0.9,font=2)
  
  write.table("Summary of Lengths of Phylum-level Annotation of Contigs",tableName,row.names=F,col.names=F,quote=F,append=T,sep="\t")
  write.table(rbind(names(evals$phylumTabNt),t(evals$phylumTabNt)),tableName,row.names=F,col.names=F,quote=F,append=T,sep="\t")
  
  evals$classTab <- sort(table(popOIcontigs$class.Kraken),decreasing=T)
  
  write.table("Class-level Annotation of Contigs",tableName,row.names=F,col.names=F,quote=F,append=T,sep="\t")
  write.table(rbind(names(evals$classTab),t(evals$classTab)),tableName,row.names=F,col.names=F,quote=F,append=T,sep="\t")
  
  classTabNt <- aggregate(popOIcontigs$length,list(popOIcontigs$class.Kraken),sum)
  rownames(classTabNt) <- classTabNt$Group.1
  evals$classTabNt <- classTabNt$x
  names(evals$classTabNt) <- classTabNt$Group.1
  evals$classTabNt <- sort(evals$classTabNt,decreasing=T)
  par(mar=c(0.2,1,2.4,1))
  plab <- names(evals$classTab)
  pie(evals$classTab,labels=plab,col=brewer.pal(9,"Set1"),clockwise=T)
  mtext("Class level annotation (Kraken)",3,1,font=2)
  mtext("contig number",3,0.2,cex=0.9)
  par(mar=c(1.1,1,1.5,1))
  plab <- sort(unique(popOIcontigs$class.Kraken))[order(classTabNt$x,decreasing=T)]
  pie(evals$classTabNt,labels=plab,col=brewer.pal(9,"Set1"),clockwise=T)
  mtext("contig length",3,0.2,cex=0.9)
  
  write.table("Summary of Lengths of Class-level Annotation of Contigs",tableName,row.names=F,col.names=F,quote=F,append=T,sep="\t")
  write.table(rbind(names(evals$classTabNt),t(evals$classTabNt)),tableName,row.names=F,col.names=F,quote=F,append=T,sep="\t")
  
  evals$orderTab <- sort(table(popOIcontigs$order.Kraken),decreasing=T)
  
  write.table("Order-level Annotation of Contigs",tableName,row.names=F,col.names=F,quote=F,append=T,sep="\t")
  write.table(rbind(names(evals$orderTab),t(evals$orderTab)),tableName,row.names=F,col.names=F,quote=F,append=T,sep="\t")
  
  orderTabNt <- aggregate(popOIcontigs$length,list(popOIcontigs$order.Kraken),sum)
  rownames(orderTabNt) <- orderTabNt$Group.1
  evals$orderTabNt <- orderTabNt$x
  names(evals$orderTabNt) <- orderTabNt$Group.1
  evals$orderTabNt <- sort(evals$orderTabNt,decreasing=T)
  par(mar=c(0.2,1,2.4,1))
  plab <- names(evals$orderTab)
  pie(evals$orderTab,labels=plab,col=brewer.pal(9,"Set1"),clockwise=T)
  mtext("Order level annotation (Kraken)",3,1,font=2)
  mtext("contig number",3,0.2,cex=0.9)
  par(mar=c(1.1,1,1.5,1))
  plab <- sort(unique(popOIcontigs$order.Kraken))[order(orderTabNt$x,decreasing=T)]
  pie(evals$orderTabNt,labels=plab,col=brewer.pal(9,"Set1"),clockwise=T)
  mtext("contig length",3,0.2,cex=0.9)
  
  write.table("Summary of Lengths of Order-level Annotation of Contigs",tableName,row.names=F,col.names=F,quote=F,append=T,sep="\t")
  write.table(rbind(names(evals$orderTabNt),t(evals$orderTabNt)),tableName,row.names=F,col.names=F,quote=F,append=T,sep="\t")
  
  evals$familyTab <- sort(table(popOIcontigs$family.Kraken),decreasing=T)
  
  write.table("Family-level Annotation of Contigs",tableName,row.names=F,col.names=F,quote=F,append=T,sep="\t")
  write.table(rbind(names(evals$familyTab),t(evals$familyTab)),tableName,row.names=F,col.names=F,quote=F,append=T,sep="\t")
  
  familyTabNt <- aggregate(popOIcontigs$length,list(popOIcontigs$family.Kraken),sum)
  rownames(familyTabNt) <- familyTabNt$Group.1
  evals$familyTabNt <- familyTabNt$x
  names(evals$familyTabNt) <- familyTabNt$Group.1
  evals$familyTabNt <- sort(evals$familyTabNt,decreasing=T)
  par(mar=c(0.2,1,2.4,1))
  plab <- names(evals$familyTab)
  pie(evals$familyTab,labels=plab,col=brewer.pal(9,"Set1"),clockwise=T)
  mtext("Family level annotation (Kraken)",3,1,font=2)
  mtext("contig number",3,0.2,cex=0.9)
  par(mar=c(1.1,1,1.5,1))
  plab <- sort(unique(popOIcontigs$family.Kraken))[order(familyTabNt$x,decreasing=T)]
  pie(evals$familyTabNt,labels=plab,col=brewer.pal(9,"Set1"),clockwise=T)
  mtext("contig length",3,0.2,cex=0.9)
  
  write.table("Summary of Lengths of Family-level Annotation of Contigs",tableName,row.names=F,col.names=F,quote=F,append=T,sep="\t")
  write.table(rbind(names(evals$familyTabNt),t(evals$familyTabNt)),tableName,row.names=F,col.names=F,quote=F,append=T,sep="\t")
  
  evals$genusTab <- sort(table(popOIcontigs$genus.Kraken),decreasing=T)
  
  write.table("Genus-level Annotation of Contigs",tableName,row.names=F,col.names=F,quote=F,append=T,sep="\t")
  write.table(rbind(names(evals$genusTab),t(evals$genusTab)),tableName,row.names=F,col.names=F,quote=F,append=T,sep="\t")
  
  genusTabNt <- aggregate(popOIcontigs$length,list(popOIcontigs$genus.Kraken),sum)
  rownames(genusTabNt) <- genusTabNt$Group.1
  evals$genusTabNt <- genusTabNt$x
  names(evals$genusTabNt) <- genusTabNt$Group.1
  evals$genusTabNt <- sort(evals$genusTabNt,decreasing=T)
  par(mar=c(0.2,1,2.4,1))
  plab <- names(evals$genusTab)
  pie(evals$genusTab,labels=plab,col=brewer.pal(9,"Set1"),clockwise=T)
  mtext("Genus level annotation (Kraken)",3,1,font=2)
  mtext("contig number",3,0.2,cex=0.9)
  par(mar=c(1.1,1,1.5,1))
  plab <- sort(unique(popOIcontigs$genus.Kraken))[order(genusTabNt$x,decreasing=T)]
  pie(evals$genusTabNt,labels=plab,col=brewer.pal(9,"Set1"),clockwise=T)
  mtext("contig length",3,0.2,cex=0.9)
  
  write.table("Summary of Lengths of Genus-level Annotation of Contigs",tableName,row.names=F,col.names=F,quote=F,append=T,sep="\t")
  write.table(rbind(names(evals$genusTabNt),t(evals$genusTabNt)),tableName,row.names=F,col.names=F,quote=F,append=T,sep="\t")
  
  evals$speciesTab <- sort(table(popOIcontigs$species.Kraken),decreasing=T)
  
  write.table("Species-level Annotation of Contigs",tableName,row.names=F,col.names=F,quote=F,append=T,sep="\t")
  write.table(rbind(names(evals$speciesTab),t(evals$speciesTab)),tableName,row.names=F,col.names=F,quote=F,append=T,sep="\t")
  
  speciesTabNt <- aggregate(popOIcontigs$length,list(popOIcontigs$species.Kraken),sum)
  rownames(speciesTabNt) <- speciesTabNt$Group.1
  evals$speciesTabNt <- speciesTabNt$x
  names(evals$speciesTabNt) <- speciesTabNt$Group.1
  evals$speciesTabNt <- sort(evals$speciesTabNt,decreasing=T)
  par(mar=c(0.2,1.5,2.4,1.5))
  plab <- names(evals$speciesTab)
  pie(evals$speciesTab,labels=plab,col=brewer.pal(9,"Set1"),clockwise=T)
  mtext("Species level annotation (Kraken)",3,1,font=2)
  mtext("contig number",3,0.2,cex=0.9)
  par(mar=c(1.1,1.5,1.5,1.5))
  plab <- sort(unique(popOIcontigs$species.Kraken))[order(speciesTabNt$x,decreasing=T)]
  pie(evals$speciesTabNt,labels=plab,col=brewer.pal(9,"Set1"),clockwise=T)
  mtext("contig length",3,0.2,cex=0.9)
  
  write.table("Summary of Lengths of Species-level Annotation of Contigs",tableName,row.names=F,col.names=F,quote=F,append=T,sep="\t")
  write.table(rbind(names(evals$speciesTabNt),t(evals$speciesTabNt)),tableName,row.names=F,col.names=F,quote=F,append=T,sep="\t")
  
  dev.off()
  return(evals)
}
makeBedContig <- function(contigTab){
  bed <- contigTab[,c("contig","length")]
  bed$start <- rep(0,nrow(bed))
  bed <- bed[,c(1,3,2)]
}
makeBedGene <- function(geneTab){
  bed <- geneTab[,c("contig","start","end","gene","sense")]
  bed$start <- bed$start-1
  bed$score <- rep(0,nrow(bed))
  bed <- bed[,c(1:4,6,5)]
}

muClus <- function(clusterName,cRes,resFile){
  dInfo <- contigInfo[cRes$class==clusterName,]
  if(dip.test(log10(1+dInfo[,4][dInfo$essentialGene!="notEssential"]))$p.value<0.05){
    cuto <- find.cutoff(log10(1+dInfo[,4][dInfo$essentialGene!="notEssential"]))
    write.table(t(c(clusterName,cuto)),resFile,append=T,row.names=F,col.names=F,quote=F,sep="\t")
	subs1 <- dInfo[log10(1+dInfo[,4])<cuto,]
    num1 <- length(unlist(sapply(subs1$essentialGene[subs1$essentialGene!="notEssential"],function(x) unlist(strsplit(x,split=";")))))
    uni1 <- length(unique(unlist(sapply(subs1$essentialGene[subs1$essentialGene!="notEssential"],function(x) unlist(strsplit(x,split=";"))))))
    if(uni1==0) {
      cRes$class[cRes$contig %in% subs1$contig] <- paste(gsub("D","E",cRes$class[cRes$contig %in% subs1$contig]),1,sep=".")
    } else if(num1/uni1<=1.2){
      cRes$class[cRes$contig %in% subs1$contig] <- paste(gsub("D","C",cRes$class[cRes$contig %in% subs1$contig]),1,sep=".")
    } else {
      cRes$class[cRes$contig %in% subs1$contig] <- paste(cRes$class[cRes$contig %in% subs1$contig],1,sep=".")
      cRes <- muClus(unique(cRes$class[cRes$contig %in% subs1$contig]),cRes,resFile)
    }
    subs2 <- dInfo[log10(1+dInfo[,4])>=cuto,]
    num2 <- length(unlist(sapply(subs2$essentialGene[subs2$essentialGene!="notEssential"],function(x) unlist(strsplit(x,split=";")))))
    uni2 <- length(unique(unlist(sapply(subs2$essentialGene[subs2$essentialGene!="notEssential"],function(x) unlist(strsplit(x,split=";"))))))
    if(uni2==0) {
      cRes$class[cRes$contig %in% subs2$contig] <- paste(gsub("D","E",cRes$class[cRes$contig %in% subs2$contig]),2,sep=".")
    } else if(num2/uni2<=1.2){
      cRes$class[cRes$contig %in% subs2$contig] <- paste(gsub("D","C",cRes$class[cRes$contig %in% subs2$contig]),2,sep=".")
    } else {
      cRes$class[cRes$contig %in% subs2$contig] <- paste(cRes$class[cRes$contig %in% subs2$contig],2,sep=".")
      cRes <- muClus(unique(cRes$class[cRes$contig %in% subs2$contig]),cRes,resFile)
    }
  } else {
    cRes$class[cRes$contig %in% dInfo$contig] <- gsub("^.","B",cRes$class[cRes$contig %in% dInfo$contig])
  }
  return(cRes)
}

coco <- ncol(contigInfo)-c(1:0)
skn <- sort(knn.dist(contigInfo[,coco],4)[,4])
sdkn <- runsd(skn,10)
est <- sort(skn)[min(which(sdkn>quantile(sdkn,0.975)&skn>median(skn)))]
write.table(t(c("scan","reachabilityDist")),"reachabilityDistanceEstimates.tsv",row.names=F,col.names=F,quote=F,sep="\t")
write.table(t(c("first",est)),"reachabilityDistanceEstimates.tsv",row.names=F,col.names=F,quote=F,sep="\t",append=T)
cdb <- dbscan(contigInfo[,coco],est,10)

plot(contigInfo[,coco],pch=16,cex=0.25,ann=F,axes=F)
j <- 1
cdbTab <- data.frame("cluster"=names(table(cdb$cluster)),"contigs"=0,"numberEss"=0,"uniqueEss"=0,stringsAsFactors=F)
for(i in names(table(cdb$cluster))) {
  points(contigInfo[cdb$cluster==i,coco],pch=16,cex=0.4,col=c("grey",colors(distinct=T))[j])
  cdbTab$contigs[cdbTab$cluster==i] <- table(cdb$cluster)[names(table(cdb$cluster))==i]
  cdbTab$numberEss[cdbTab$cluster==i] <- length(unlist(sapply(contigInfo$essentialGene[contigInfo$essentialGene!="notEssential"&cdb$cluster==i],function(x) unlist(strsplit(x,split=";")))))
  cdbTab$uniqueEss[cdbTab$cluster==i] <- length(unique(unlist(sapply(contigInfo$essentialGene[contigInfo$essentialGene!="notEssential"&cdb$cluster==i],function(x) unlist(strsplit(x,split=";"))))))
  j<-j+1
}
box()
write.table(cdbTab,"clusterFiles/clusterFirstScan.tsv",sep="\t",row.names=F,quote=F)
write.table(t(c("clusterName","cutoff")),"bimodalClusterCutoffs.tsv",sep="\t",row.names=F,col.names=F,quote=F)
clusterRes <- data.frame("contig"=contigInfo$contig,"class"="x",stringsAsFactors=F)
clusterRes$class[cdb$cluster==0] <- "N"

emptyClus <- cdbTab$cluster[cdbTab$uniqueEss==0&cdbTab$cluster!=0]
clusterRes$class[cdb$cluster %in% emptyClus] <- paste("E",cdb$cluster[cdb$cluster %in% emptyClus],sep="")

closeClus <- cdbTab$cluster[cdbTab$numberEss/cdbTab$uniqueEss<=1.2&cdbTab$cluster!=0]
clusterRes$class[cdb$cluster %in% closeClus] <- paste("C",cdb$cluster[cdb$cluster %in% closeClus],sep="")

duClus <- cdbTab$cluster[cdbTab$numberEss/cdbTab$uniqueEss>1.2&cdbTab$cluster!=0]
clusterRes$class[cdb$cluster %in% duClus] <- paste("D",cdb$cluster[cdb$cluster %in% duClus],sep="")

write.table(t(c("cluster","cutoff")),"bimodalClusterCutoffs.tsv",sep="\t",col.names=F,row.names=F,quote=F)
for(ds in unique(clusterRes$class[grep("D",clusterRes$class)])){
	clusterRes <- muClus(ds,clusterRes,"bimodalClusterCutoffs.tsv")
}

for(bb in unique(clusterRes$class[grep("B",clusterRes$class)])){
	bbInfo <- contigInfo[clusterRes$class==bb,]
  sknBB <- sort(knn.dist(bbInfo[,coco],4)[,4])
  if(length(sknBB>10)){
    sdknBB <- runsd(sknBB,10)
    estBB <- sort(sknBB)[min(which(sdknBB>quantile(sdknBB,0.975)&sknBB>median(sknBB)))]
    if(is.na(estBB)) estBB <- sort(sknBB)[max(sdknBB)]

  write.table(t(c(bb,estBB)),"reachabilityDistanceEstimates.tsv",row.names=F,col.names=F,quote=F,sep="\t",append=T)
  BBcdb <- dbscan(bbInfo[,coco],estBB,15)
  
  plot(bbInfo[,coco],pch=16,cex=0.25,ann=F)
  j <- 1
BBcdbTab <- data.frame("cluster"=names(table(BBcdb$cluster)),"contigs"=0,"numberEss"=0,"uniqueEss"=0,stringsAsFactors=F)
    for(i in names(table(BBcdb$cluster))) {
      points(bbInfo[BBcdb$cluster==i,coco],pch=16,cex=0.4,col=c("grey",colors(distinct=T))[j])
      BBcdbTab$contigs[BBcdbTab$cluster==i] <- table(BBcdb$cluster)[names(table(BBcdb$cluster))==i]
      BBcdbTab$numberEss[BBcdbTab$cluster==i] <- length(unlist(sapply(bbInfo$essentialGene[bbInfo$essentialGene!="notEssential"&BBcdb$cluster==i],function(x) unlist(strsplit(x,split=";")))))
      BBcdbTab$uniqueEss[BBcdbTab$cluster==i] <- length(unique(unlist(sapply(bbInfo$essentialGene[bbInfo$essentialGene!="notEssential"&BBcdb$cluster==i],function(x) unlist(strsplit(x,split=";"))))))
      j<-j+1
    }
    write.table(BBcdbTab,paste("clusterFiles/blob1Cluster",bb,".tsv",sep=""),sep="\t",row.names=F,quote=F)
    
    BBclusterRes <- data.frame("contig"=bbInfo$contig,"class"="x",stringsAsFactors=F)
    BBclusterRes$class[BBcdb$cluster==0] <- "N"
    clusterRes$class[clusterRes$contig %in% BBclusterRes$contig[BBclusterRes$class=="N"]] <- "N"
    
    BBemptyClus <- BBcdbTab$cluster[BBcdbTab$uniqueEss==0&BBcdbTab$cluster!=0]
    BBclusterRes$class[BBcdb$cluster %in% BBemptyClus] <- paste("E",BBcdb$cluster[BBcdb$cluster %in% BBemptyClus],sep="")
    BBemptyCont <- BBclusterRes$contig[BBcdb$cluster %in% BBemptyClus]
    clusterRes$class[clusterRes$contig %in% BBemptyCont] <- paste(gsub("B","E",bb),
                                                                  gsub("E","",unlist(sapply(clusterRes$contig[clusterRes$contig %in% BBemptyCont],function(x)BBclusterRes$class[BBclusterRes$contig==x]))),
                                                                  sep=".")
    
    BBcloseClus <- BBcdbTab$cluster[BBcdbTab$numberEss/BBcdbTab$uniqueEss<=1.2&BBcdbTab$cluster!=0]
    BBclusterRes$class[BBcdb$cluster %in% BBcloseClus] <- paste("C",BBcdb$cluster[BBcdb$cluster %in% BBcloseClus],sep="")
    BBcloseCont <- BBclusterRes$contig[BBcdb$cluster %in% BBcloseClus]
    clusterRes$class[clusterRes$contig %in% BBcloseCont] <- paste(gsub("B","C",bb),
                                                                  gsub("C","",unlist(sapply(clusterRes$contig[clusterRes$contig %in% BBcloseCont],function(x)BBclusterRes$class[BBclusterRes$contig==x]))),
                                                                  sep=".")
    
    BBduClus <- BBcdbTab$cluster[BBcdbTab$numberEss/BBcdbTab$uniqueEss>1.2&BBcdbTab$cluster!=0]
    BBclusterRes$class[BBcdb$cluster %in% BBduClus] <- paste("D",BBcdb$cluster[BBcdb$cluster %in% BBduClus],sep="")
    BBduCont <- BBclusterRes$contig[BBcdb$cluster %in% BBduClus]
    uniD <- unique(paste(gsub("B","D",bb),
                         gsub("D","",unlist(sapply(clusterRes$contig[clusterRes$contig %in% BBduCont],function(x)BBclusterRes$class[BBclusterRes$contig==x]))),
                         sep="."))
    clusterRes$class[clusterRes$contig %in% BBduCont] <- paste(gsub("B","D",bb),
                                                               gsub("D","",unlist(sapply(clusterRes$contig[clusterRes$contig %in% BBduCont],function(x)BBclusterRes$class[BBclusterRes$contig==x]))),
                                                               sep=".")
    for(bds in uniD){
      clusterRes <- muClus(bds,clusterRes,"bimodalClusterCutoffs.tsv")
    }
  }
}

for(bb in unique(clusterRes$class[grep("B",clusterRes$class)])){
        bbInfo <- contigInfo[clusterRes$class==bb,]
  sknBB <- sort(knn.dist(bbInfo[,coco],4)[,4])
  if(length(sknBB>10)){
    sdknBB <- runsd(sknBB,10)
    estBB <- sort(sknBB)[min(which(sdknBB>quantile(sdknBB,0.975)&sknBB>median(sknBB)))]
    if(is.na(estBB)) estBB <- sort(sknBB)[max(sdknBB)]

  write.table(t(c(bb,estBB)),"reachabilityDistanceEstimates.tsv",row.names=F,col.names=F,quote=F,sep="\t",append=T)
  BBcdb <- dbscan(bbInfo[,coco],estBB,12)

  plot(bbInfo[,coco],pch=16,cex=0.25,ann=F)
  j <- 1
BBcdbTab <- data.frame("cluster"=names(table(BBcdb$cluster)),"contigs"=0,"numberEss"=0,"uniqueEss"=0,stringsAsFactors=F)
    for(i in names(table(BBcdb$cluster))) {
      points(bbInfo[BBcdb$cluster==i,coco],pch=16,cex=0.4,col=c("grey",colors(distinct=T))[j])
      BBcdbTab$contigs[BBcdbTab$cluster==i] <- table(BBcdb$cluster)[names(table(BBcdb$cluster))==i]
      BBcdbTab$numberEss[BBcdbTab$cluster==i] <- length(unlist(sapply(bbInfo$essentialGene[bbInfo$essentialGene!="notEssential"&BBcdb$cluster==i],function(x) unlist(strsplit(x,split=";")))))
      BBcdbTab$uniqueEss[BBcdbTab$cluster==i] <- length(unique(unlist(sapply(bbInfo$essentialGene[bbInfo$essentialGene!="notEssential"&BBcdb$cluster==i],function(x) unlist(strsplit(x,split=";"))))))
      j<-j+1
    }
    write.table(BBcdbTab,paste("clusterFiles/blob2Cluster",bb,".tsv",sep=""),sep="\t",row.names=F,quote=F)

    BBclusterRes <- data.frame("contig"=bbInfo$contig,"class"="x",stringsAsFactors=F)
    BBclusterRes$class[BBcdb$cluster==0] <- "N"
    clusterRes$class[clusterRes$contig %in% BBclusterRes$contig[BBclusterRes$class=="N"]] <- "N"
    BBemptyClus <- BBcdbTab$cluster[BBcdbTab$uniqueEss==0&BBcdbTab$cluster!=0]
    BBclusterRes$class[BBcdb$cluster %in% BBemptyClus] <- paste("E",BBcdb$cluster[BBcdb$cluster %in% BBemptyClus],sep="")
    BBemptyCont <- BBclusterRes$contig[BBcdb$cluster %in% BBemptyClus]
    clusterRes$class[clusterRes$contig %in% BBemptyCont] <- paste(gsub("B","E",bb),
                                                                  gsub("E","",unlist(sapply(clusterRes$contig[clusterRes$contig %in% BBemptyCont],function(x)BBclusterRes$class[BBclusterRes$contig==x]))),
                                                                  sep=".")

    BBcloseClus <- BBcdbTab$cluster[BBcdbTab$numberEss/BBcdbTab$uniqueEss<=1.2&BBcdbTab$cluster!=0]
    BBclusterRes$class[BBcdb$cluster %in% BBcloseClus] <- paste("C",BBcdb$cluster[BBcdb$cluster %in% BBcloseClus],sep="")
    BBcloseCont <- BBclusterRes$contig[BBcdb$cluster %in% BBcloseClus]
    clusterRes$class[clusterRes$contig %in% BBcloseCont] <- paste(gsub("B","C",bb),
                                                                  gsub("C","",unlist(sapply(clusterRes$contig[clusterRes$contig %in% BBcloseCont],function(x)BBclusterRes$class[BBclusterRes$contig==x]))),
                                                                  sep=".")

    BBduClus <- BBcdbTab$cluster[BBcdbTab$numberEss/BBcdbTab$uniqueEss>1.2&BBcdbTab$cluster!=0]
    BBclusterRes$class[BBcdb$cluster %in% BBduClus] <- paste("D",BBcdb$cluster[BBcdb$cluster %in% BBduClus],sep="")
    BBduCont <- BBclusterRes$contig[BBcdb$cluster %in% BBduClus]
    uniD <- unique(paste(gsub("B","D",bb),
                         gsub("D","",unlist(sapply(clusterRes$contig[clusterRes$contig %in% BBduCont],function(x)BBclusterRes$class[BBclusterRes$contig==x]))),
                         sep="."))
    clusterRes$class[clusterRes$contig %in% BBduCont] <- paste(gsub("B","D",bb),
                                                               gsub("D","",unlist(sapply(clusterRes$contig[clusterRes$contig %in% BBduCont],function(x)BBclusterRes$class[BBclusterRes$contig==x]))),
                                                               sep=".")
    for(bds in uniD){
      clusterRes <- muClus(bds,clusterRes,"bimodalClusterCutoffs.tsv")
    }
  }
}

for(bb in unique(clusterRes$class[grep("B",clusterRes$class)])){
        bbInfo <- contigInfo[clusterRes$class==bb,]
  sknBB <- sort(knn.dist(bbInfo[,coco],4)[,4])
  if(length(sknBB>10)){
    sdknBB <- runsd(sknBB,10)
    estBB <- sort(sknBB)[min(which(sdknBB>quantile(sdknBB,0.975)&sknBB>median(sknBB)))]
    if(is.na(estBB)) estBB <- sort(sknBB)[max(sdknBB)]

  write.table(t(c(bb,estBB)),"reachabilityDistanceEstimates.tsv",row.names=F,col.names=F,quote=F,sep="\t",append=T)
  BBcdb <- dbscan(bbInfo[,coco],estBB,15)

  plot(bbInfo[,coco],pch=16,cex=0.25,ann=F)
  j <- 1
BBcdbTab <- data.frame("cluster"=names(table(BBcdb$cluster)),"contigs"=0,"numberEss"=0,"uniqueEss"=0,stringsAsFactors=F)
    for(i in names(table(BBcdb$cluster))) {
      points(bbInfo[BBcdb$cluster==i,coco],pch=16,cex=0.4,col=c("grey",colors(distinct=T))[j])
      BBcdbTab$contigs[BBcdbTab$cluster==i] <- table(BBcdb$cluster)[names(table(BBcdb$cluster))==i]
      BBcdbTab$numberEss[BBcdbTab$cluster==i] <- length(unlist(sapply(bbInfo$essentialGene[bbInfo$essentialGene!="notEssential"&BBcdb$cluster==i],function(x) unlist(strsplit(x,split=";")))))
      BBcdbTab$uniqueEss[BBcdbTab$cluster==i] <- length(unique(unlist(sapply(bbInfo$essentialGene[bbInfo$essentialGene!="notEssential"&BBcdb$cluster==i],function(x) unlist(strsplit(x,split=";"))))))
      j<-j+1
    }
    write.table(BBcdbTab,paste("clusterFiles/blob3Cluster",bb,".tsv",sep=""),sep="\t",row.names=F,quote=F)

    BBclusterRes <- data.frame("contig"=bbInfo$contig,"class"="x",stringsAsFactors=F)
    BBclusterRes$class[BBcdb$cluster==0] <- "N"
    clusterRes$class[clusterRes$contig %in% BBclusterRes$contig[BBclusterRes$class=="N"]] <- "N"
    BBemptyClus <- BBcdbTab$cluster[BBcdbTab$uniqueEss==0&BBcdbTab$cluster!=0]
    BBclusterRes$class[BBcdb$cluster %in% BBemptyClus] <- paste("E",BBcdb$cluster[BBcdb$cluster %in% BBemptyClus],sep="")
    BBemptyCont <- BBclusterRes$contig[BBcdb$cluster %in% BBemptyClus]
    clusterRes$class[clusterRes$contig %in% BBemptyCont] <- paste(gsub("B","E",bb),
                                                                  gsub("E","",unlist(sapply(clusterRes$contig[clusterRes$contig %in% BBemptyCont],function(x)BBclusterRes$class[BBclusterRes$contig==x]))),
                                                                  sep=".")

    BBcloseClus <- BBcdbTab$cluster[BBcdbTab$numberEss/BBcdbTab$uniqueEss<=1.2&BBcdbTab$cluster!=0]
    BBclusterRes$class[BBcdb$cluster %in% BBcloseClus] <- paste("C",BBcdb$cluster[BBcdb$cluster %in% BBcloseClus],sep="")
    BBcloseCont <- BBclusterRes$contig[BBcdb$cluster %in% BBcloseClus]
    clusterRes$class[clusterRes$contig %in% BBcloseCont] <- paste(gsub("B","C",bb),
                                                                  gsub("C","",unlist(sapply(clusterRes$contig[clusterRes$contig %in% BBcloseCont],function(x)BBclusterRes$class[BBclusterRes$contig==x]))),
                                                                  sep=".")

    BBduClus <- BBcdbTab$cluster[BBcdbTab$numberEss/BBcdbTab$uniqueEss>1.2&BBcdbTab$cluster!=0]
    BBclusterRes$class[BBcdb$cluster %in% BBduClus] <- paste("D",BBcdb$cluster[BBcdb$cluster %in% BBduClus],sep="")
    BBduCont <- BBclusterRes$contig[BBcdb$cluster %in% BBduClus]
    uniD <- unique(paste(gsub("B","D",bb),
                         gsub("D","",unlist(sapply(clusterRes$contig[clusterRes$contig %in% BBduCont],function(x)BBclusterRes$class[BBclusterRes$contig==x]))),
                         sep="."))
    clusterRes$class[clusterRes$contig %in% BBduCont] <- paste(gsub("B","D",bb),
                                                               gsub("D","",unlist(sapply(clusterRes$contig[clusterRes$contig %in% BBduCont],function(x)BBclusterRes$class[BBclusterRes$contig==x]))),
                                                               sep=".")
    for(bds in uniD){
      clusterRes <- muClus(bds,clusterRes,"bimodalClusterCutoffs.tsv")
    }
  }
}

for(clus in grep("C",unique(clusterRes$class),value=T)){
  uni <- length(unique(geneInfo$essentialGene[geneInfo$essentialGene!="notEssential"&geneInfo$contig %in% contigInfo$contig[clusterRes$class==clus]]))
  ess <- length(geneInfo$essentialGene[geneInfo$essentialGene!="notEssential"&geneInfo$contig %in% contigInfo$contig[clusterRes$class==clus]])
  if(uni>100&ess<115){
    clusterRes$class[clusterRes$class == clus] <- gsub("C","P",clus)
  } else if(uni>71){
    clusterRes$class[clusterRes$class == clus] <- gsub("C","G",clus)
  } else if(uni>51){
    clusterRes$class[clusterRes$class == clus] <- gsub("C","O",clus)
  } else if(uni>31){
    clusterRes$class[clusterRes$class == clus] <- gsub("C","L",clus)
  }
}

write.table(clusterRes,"contigs2clusters.tsv",sep="\t",row.names=F,quote=F)
essPal <- colorRampPalette(brewer.pal(11,"Spectral"))(111)[111:1]
plot(contigInfo[,coco],pch=16,cex=0.25,ann=F,axes=F,bty="o")
clusterInfo <- data.frame("cluster"=unique(clusterRes$class),"length"=0,"contigs"=0,"aveCov"=0,"varPerMB"=0,"genes"=0,"completeGenes"=0,"expressedGenes"=0,
                          "uniqueEss"=0,"numEss"=0,"catsKO"=0,"catsMetaCyc"=0,"catsSwissProt"=0,"catsPfam"=0,"catsTIGR"=0,
                          "annotated"=0,"annotatedKO"=0,"annotatedMetaCyc"=0,"annotatedSwissProt"=0,"annotatedPfam"=0,"annotatedTIGR"=0,
                          "phylaKraken"="","generaKraken"="","phylaAmphora"="","generaAmphora"="","phylaMotuPresent"="",
                          "generaMotuPresent"="",stringsAsFactors=F)
for(i in 1:nrow(clusterInfo)){
  clus <- clusterInfo$cluster[i]
  clusterInfo$length[i] <- sum(contigInfo$length[contigInfo$contig %in% clusterRes$contig[clusterRes$class==clus]])
  clusterInfo$contigs[i] <- nrow(contigInfo[contigInfo$contig %in% clusterRes$contig[clusterRes$class==clus],])
  clusterInfo$aveCov[i] <- weighted.mean(contigInfo[contigInfo$contig %in% clusterRes$contig[clusterRes$class==clus],4],
                                         contigInfo$length[contigInfo$contig %in% clusterRes$contig[clusterRes$class==clus]])
  clusterInfo$varPerMB[i] <- weighted.mean(contigInfo$varPerMB[contigInfo$contig %in% clusterRes$contig[clusterRes$class==clus]],
                                           contigInfo$length[contigInfo$contig %in% clusterRes$contig[clusterRes$class==clus]])
  clusterInfo$genes[i] <- nrow(geneInfo[geneInfo$contig %in% clusterRes$contig[clusterRes$class==clus],])
  clusterInfo$completeGenes[i] <- nrow(geneInfo[geneInfo$contig %in% clusterRes$contig[clusterRes$class==clus]&
                                                  geneInfo$completeness=="complete",])
  clusterInfo$expressedGenes[i] <- nrow(geneInfo[geneInfo$contig %in% clusterRes$contig[clusterRes$class==clus]&
                                                   geneInfo$aveCovRNAfw>=0.01,])
  clusterInfo$uniqueEss[i] <- length(unique(geneInfo$essentialGene[geneInfo$contig %in% clusterRes$contig[clusterRes$class==clus]&
                                                                     geneInfo$essentialGene!="notEssential"]))
  clusterInfo$numEss[i] <- length(geneInfo$essentialGene[geneInfo$contig %in% clusterRes$contig[clusterRes$class==clus]&
                                                           geneInfo$essentialGene!="notEssential"])
  if(clus=="N") pc <- "grey" else if(grepl("E",clus)) pc <- "black" else if(grepl("B",clus)) pc <- brewer.pal(8,"Set3")[8] else{
    pc <- essPal[clusterInfo$uniqueEss[i]]
  }
  points(contigInfo[clusterRes$class==clus,coco],cex=0.3,pch=16,col=pc)
  clusterInfo$catsKO[i] <- length(unique(unlist(sapply(geneInfo$KO[geneInfo$contig %in% clusterRes$contig[clusterRes$class==clus]&
                                                                     geneInfo$KO!="unknown"],
                                                       function(x) unlist(strsplit(x,split=";"))))))
  clusterInfo$catsMetaCyc[i] <- length(unique(unlist(sapply(geneInfo$metaCycID[geneInfo$contig %in% clusterRes$contig[clusterRes$class==clus]&
                                                                                 geneInfo$metaCycID!="unknown"],
                                                            function(x) unlist(strsplit(x,split=";"))))))
  clusterInfo$catsSwissProt[i] <- length(unique(unlist(sapply(geneInfo$swissprotEC[geneInfo$contig %in% clusterRes$contig[clusterRes$class==clus]&
                                                                                     geneInfo$swissprotEC!="unknown"],
                                                              function(x) unlist(strsplit(x,split=";"))))))
  clusterInfo$catsPfam[i] <- length(unique(unlist(sapply(geneInfo$pfamID[geneInfo$contig %in% clusterRes$contig[clusterRes$class==clus]&
                                                                           geneInfo$pfamID!="unknown"],
                                                         function(x) unlist(strsplit(x,split=";"))))))
  clusterInfo$catsTIGR[i] <- length(unique(unlist(sapply(geneInfo$tigrID[geneInfo$contig %in% clusterRes$contig[clusterRes$class==clus]&
                                                                           geneInfo$tigrID!="unknown"],
                                                         function(x) unlist(strsplit(x,split=";"))))))
  clusterInfo$annotated[i] <- nrow(geneInfo[geneInfo$contig %in% clusterRes$contig[clusterRes$class==clus]&
                                              (geneInfo$KO!="unknown"|geneInfo$metaCycID!="unknown"|geneInfo$swissprotEC!="unknown"|
                                              geneInfo$pfamID!="unknown"|geneInfo$tigrID!="unknown"),])
  clusterInfo$annotatedKO[i] <- nrow(geneInfo[geneInfo$contig %in% clusterRes$contig[clusterRes$class==clus]&
                                                geneInfo$KO!="unknown",])
  clusterInfo$annotatedMetaCyc[i] <- nrow(geneInfo[geneInfo$contig %in% clusterRes$contig[clusterRes$class==clus]&
                                                     geneInfo$metaCycID!="unknown",])
  clusterInfo$annotatedSwissProt[i] <- nrow(geneInfo[geneInfo$contig %in% clusterRes$contig[clusterRes$class==clus]&
                                                       geneInfo$swissprotEC!="unknown",])
  clusterInfo$annotatedPfam[i] <- nrow(geneInfo[geneInfo$contig %in% clusterRes$contig[clusterRes$class==clus]&
                                                  geneInfo$pfamID!="unknown",])
  clusterInfo$annotatedTIGR[i] <- nrow(geneInfo[geneInfo$contig %in% clusterRes$contig[clusterRes$class==clus]&
                                                  geneInfo$tigrID!="unknown",])
  taxt <- sort(table(contigInfo$phylum.Kraken[contigInfo$contig %in% clusterRes$contig[clusterRes$class==clus]]),decreasing=T)
  clusterInfo$phylaKraken[i] <- paste(paste(names(taxt),"(",taxt,")",sep=""),sep=";",collapse=";")
  rm(taxt)
  taxt <- sort(table(contigInfo$genus.Kraken[contigInfo$contig %in% clusterRes$contig[clusterRes$class==clus]]),decreasing=T)
  clusterInfo$generaKraken[i] <- paste(paste(names(taxt),"(",taxt,")",sep=""),sep=";",collapse=";")
  rm(taxt)
  taxt <- sort(table(unlist(sapply(geneInfo$phylum.Amphora[geneInfo$contig %in% clusterRes$contig[clusterRes$class==clus]&geneInfo$marker.Amphora!="unknown"],function(x)unlist(strsplit(x,split="\\("))[1]))),decreasing=T)
  clusterInfo$phylaAmphora[i] <- paste(paste(names(taxt),"(",taxt,")",sep=""),sep=";",collapse=";")
  rm(taxt)
  taxt <- sort(table(unlist(sapply(geneInfo$genus.Amphora[geneInfo$contig %in% clusterRes$contig[clusterRes$class==clus]&geneInfo$marker.Amphora!="unknown"],function(x)unlist(strsplit(x,split="\\("))[1]))),decreasing=T)
  clusterInfo$generaAmphora[i] <- paste(paste(names(taxt),"(",taxt,")",sep=""),sep=";",collapse=";")
  rm(taxt)
  taxt <- sort(table(geneInfo$Phylum.mOTUpresent[geneInfo$contig %in% clusterRes$contig[clusterRes$class==clus]&geneInfo$MarkerGene.mOTUpresent!="unknown"]),decreasing=T)
  clusterInfo$phylaMotuPresent[i] <- paste(paste(names(taxt),"(",taxt,")",sep=""),sep=";",collapse=";")
  rm(taxt)
  taxt <- sort(table(geneInfo$Genus.mOTUpresent[geneInfo$contig %in% clusterRes$contig[clusterRes$class==clus]&geneInfo$MarkerGene.mOTUpresent!="unknown"]),decreasing=T)
  clusterInfo$generaMotuPresent[i] <- paste(paste(names(taxt),"(",taxt,")",sep=""),sep=";",collapse=";")
  rm(taxt)
}
write.table(clusterInfo,"clusterStats.tsv",sep="\t",row.names=F,quote=F)

save.image("test.Rdata")

popInf <- list()
for(pop in c(unique(grep("P",clusterRes$class,value=T)),unique(grep("G",clusterRes$class,value=T)))) {
  popInf[[pop]] <- evalSubPop(contigInfo[contigInfo$contig %in% clusterRes$contig[clusterRes$class==pop],],geneInfo[geneInfo$contig %in% clusterRes$contig[clusterRes$class==pop],],pop,contigInfo,geneInfo)
  write.table(makeBedGene(geneInfo[geneInfo$contig %in% clusterRes$contig[clusterRes$class==pop],]),paste(LIB,"_population_",pop,".genes.bed",sep=""),quote=F,col.names=F,row.names=F,sep="\t")
  write.table(makeBedContig(contigInfo[contigInfo$contig %in% clusterRes$contig[clusterRes$class==pop],]),paste(LIB,"_population_",pop,".contigs.bed",sep=""),quote=F,col.names=F,row.names=F,sep="\t")
}

geneInf <- merge(geneInf,clusterRes,by.x=1,by.y=1,all.x=T)
geneInf[,ncol(geneInf)-2:0] <- geneInf[,ncol(geneInf)-c(0,2,1)]
colnames(geneInf)[ncol(geneInf)-2:0] <- colnames(geneInf)[ncol(geneInf)-c(0,2,1)]
geneInf$class[is.na(geneInf$class)] <- "S"
contigInf <- merge(contigInf,clusterRes,by.x=1,by.y=1,all.x=T)
contigInf[,ncol(contigInf)-2:0] <- contigInf[,ncol(contigInf)-c(0,2,1)]
colnames(contigInf)[ncol(contigInf)-2:0] <- colnames(contigInf)[ncol(contigInf)-c(0,2,1)]
contigInf$class[is.na(contigInf$class)] <- "S"

contigInfo <- contigInf[!is.na(contigInf$x),]
geneInfo <- geneInf[!is.na(geneInf$x),]

contigInfoC <- contigInfo[c(grep("P",contigInfo$class),grep("G",contigInfo$class)),]
geneInfoC <- geneInfo[c(grep("P",geneInfo$class),grep("G",geneInfo$class)),]

rm(list=c("BBbbClus","BBbbCont","BBduCont","BBcdb","BBcdbTab","BBcloseClus","BBcloseCont","BBclusterRes","BBduClus","BBemptyClus","BBemptyCont","BBmassClus","BBmassCont","bb","bbClus","bbInfo","bds","cdb","cdbTab","closeClus","clus","clusterRes","coco","dInfo","ds","duClus","emptyClus","est","estBB","ess","essPal","evalSubPop","find.cutoff","i","j","makeBedContig","makeBedGene","massClus","muClus","pc","perfClus","pop","sdkn","sdknBB","skn","sknBB","uni","uniD"))
save.image("WSclusterAll.Rdata")

rm(list=c("geneInf","contigInf"))
save.image("WSclusterLong.Rdata")

rm(list=c("geneInfo","contigInfo"))
save.image("WSclusterCompletish.Rdata")








