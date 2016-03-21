#This script contains two functions (zoomSample, zoomTaxa) which can be used to colour phylogenetic trees based on the samples of origin 
# of the phylogenetic marker genes forming its leaves. They also pull different sub-groups from the tree, 
# based on earlier taxonomic annotations or based on the samples of origin.

# Examples for usage of these functions are shown in lines 87ff.

#The script would need phylogenetic trees in Newick format (see line 98) and additional information on the chosen colours, samples of origin
# and the contig clusters the genes building the tree belong to (see line 88).

#The script is based on the file structure and sample IDs of the MuSt project and should not be expected to run on other data as is!!!

#written by Anna Heintzbuschart, August 2015

library(ape)
library(geiger)
library(RColorBrewer)

zoomSample <- function (phy, focus, subtree = FALSE, col.sub = Mcol[visFac], col.tip = compCol, markerGene = mark, ...) {
  if (!is.list(focus)) 
    focus <- list(focus)
  n <- length(focus)
  for (i in 1:n) if (is.character(focus[[i]])) 
    focus[[i]] <- which(phy$tip.label %in% focus[[i]])
  if (is.function(col.sub)) {
    col.sub <- if (deparse(substitute(col.sub)) == "grey") 
      grey(1:n/n)
    else col.sub(n)
  }
  ext <- vector("list", n)
  extcol <- vector("list", n)
  extname <- names(focus)
  for (i in 1:n) {
    ext[[i]] <- drop.tip(phy, phy$tip.label[-focus[[i]]], subtree = subtree, rooted = TRUE)
    extcol[[i]] <- unlist(sapply(ext[[i]]$tip.label,function(x) col.tip[which(phy$tip.label==x)]))
    ext[[i]]$tip.label <- gsub("M.+_","",gsub("_Contig.+","",ext[[i]]$tip.label))
  }
  nc <- round(sqrt(n)) + 1
  nr <- ceiling(sqrt(n))
  M <- matrix(0, nr, nc)
  x <- c(rep(1, nr), 2:(n + 1))
  M[1:length(x)] <- x
  layout(M, c(1, rep(3/(nc - 1), nc - 1)))
  phy$tip.label <- rep("", length(phy$tip.label))
  colo <- rep("black", dim(phy$edge)[1])
  for (i in 1:n) colo[which.edge(phy, focus[[i]])] <- col.sub[i]
  plot.phylo(phy, edge.color = colo, ...)
  mtext(markerGene,3,1)
  for (i in 1:n){
    plot.phylo(ext[[i]], edge.color = col.sub[i], tip.color = extcol[[i]])
    mtext(extname[i],3,0,col=col.sub[i],font=2)
  }
}

zoomTaxa <- function (phy, focus, subtree = FALSE, col.sub = rainbow, col.tip = sampCol, ...) {
  if (!is.list(focus)) 
    focus <- list(focus)
  n <- length(focus)
  for (i in 1:n) if (is.character(focus[[i]])) 
    focus[[i]] <- which(phy$tip.label %in% focus[[i]])
  if (is.function(col.sub)) {
    col.sub <- if (deparse(substitute(col.sub)) == "grey") 
      grey(1:n/n)
    else col.sub(n)
  }
  ext <- vector("list", n)
  extcol <- vector("list", n)
  for (i in 1:n) {
    ext[[i]] <- drop.tip(phy, phy$tip.label[-focus[[i]]], subtree = subtree, rooted = TRUE)
    extcol[[i]] <- unlist(sapply(ext[[i]]$tip.label,function(x) col.tip[which(phy$tip.label==x)]))
    ext[[i]]$tip.label <- gsub("_Contig.+","",ext[[i]]$tip.label)
  }
  nc <- round(sqrt(n)) + 1
  nr <- ceiling(sqrt(n))
  M <- matrix(0, nr, nc)
  x <- c(rep(1, nr), 2:(n + 1))
  M[1:length(x)] <- x
  layout(M, c(1, rep(3/(nc - 1), nc - 1)))
  phy$tip.label <- rep("", length(phy$tip.label))
  colo <- rep("black", dim(phy$edge)[1])
  for (i in 1:n) colo[which.edge(phy, focus[[i]])] <- col.sub[i]
  plot.phylo(phy, edge.color = colo, ...)
  for (i in 1:n){
    plot.phylo(ext[[i]], edge.color = col.sub[i], tip.color = extcol[[i]])
  }
}

#####
cI <- readRDS("clusterInf.RDS")
cis <- cI[cI$cluster!="S",]
McolA <- c(brewer.pal(7,"Blues")[3:7],brewer.pal(7,"Oranges")[3:7],brewer.pal(6,"Purples")[3:6],brewer.pal(8,"Greens")[3:8]) #all 4 families
Mcol <- McolA[-c(5,11)]
visFac <- c(rep(1:3,each=3),rep(4,times=2),rep(5,times=3),rep(6,times=2),rep(7:8,each=3),rep(9,times=2),10:13,rep(14,times=2),15,rep(16,times=2),17,rep(18,times=2))
samples <- unique(cI$sample)
compCol <- colorRampPalette(brewer.pal(11,"Spectral"))(109)[109:1]

markers <- c("rpoB","COG0012","COG0016","COG0018","COG0172","COG0215","COG0495","COG0525","COG0533","COG0541","COG0552")
for(mark in markers){
  tall <- read.tree(paste(mark,".allCompleteGenes.faa.final_tree.nw",sep=""))
  ci <- cis[cis$cluster!="N"&apply(cis[,1:2],1,function(x) paste(x,sep="_",collapse="_")) %in% 
              gsub("_Contig.+","",tall$tip.label),]

  tiplist <- list()
  for(sample in samples){
    tiplist[[sample]] <- which(gsub("_.+","",tall$tip.label)==sample)
  }
  comptip <- unlist(sapply(gsub("_Contig.+","",tall$tip.label),function(x){
    if(grepl("S",x)) "grey75" else if(grepl("N",x)) "grey50" else {
      compCol[ci$uniqueEss[apply(ci[,1:2],1,function(x) paste(x,sep="_",collapse="_")) == x]+1]
    }
  }))
  samptip <- unlist(sapply(gsub("_.+","",tall$tip.label),function(x) Mcol[visFac][which(samples==x)]))
  
  pdf(paste(mark,"perSampleTrees.pdf",sep="."),width=15,height=18,pointsize=8)
  zoomSample(tall,tiplist,col.sub=Mcol[visFac],col.tip=comptip)
  dev.off()
  
  pdf(paste(mark,"unclassPhylaTree.pdf",sep=""),width=11,height=8,pointsize=8)
  zoomTaxa(tall,which(gsub("_Contig.+","",tall$tip.label) %in% apply(ci[grep("NA",ci$phylaMotuPresent),1:2],1,
                                                                 function(x) paste(x,sep="_",collapse="_"))),
           col.tip=samptip)
  mtext("Unclassified phyla",1,1,font=2)
  mtext(mark,3,1)
  dev.off()

}


