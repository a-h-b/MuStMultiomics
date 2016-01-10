# the following functions can be used to plot modules of reconstructed networks displaying the contribution of different populations to each node
### further information, such as fold change can also be indicated
#
# needs: - the network as produced by "140630_MUST_NW.R"
#        - the result of a differential analysis of the nodes from DESeq
#        - the nodes in the module, such as the output of the BioNet functions in "140630_MUST_NW.R"
#        - a data.frame with the information on each reconstructed population-level genomes
#        - the values for every population contributing to the nodes, ie the outputs of "getNWexprMongoAllSamples.R"
#
# workflow: assuming all necessary data is in the workspace, plotModSuper() can be run
#           annaMustPrepModule() prepares the module data and creates a layout for the network
#           as a second step, it will call annaMustPlotIndiModule() which will create sample-wise plots
#           plotModSuper() takes the following arguments:
#               - "nodes" - the nodes in the topscoring module,
#               - "diffRes" - a file with the results of DESeq
#               - "indiv" - a vector with the names of the samples for which a plot should be produced
#               - "outname" - a name to be used for the output
#               - "groups" - either a vector with the group membership for each of the samples, given as 1 or 2; or just 1, if no difference should be made (defaults to 1)
#               - "ome" - a vector with "G", "T" or "P" or combinations of the letters, indicating which omic level(s) should be plotted (defaults to "T")
#               - "silent" - logical, decides whether plotting to stdout is suppressed (defaults to T, ie. only plots in .pdf are created)
#               - "compc" - a vector with colors for each level of completeness (defaults to the variable name compCol)
#               - "clusterInf" - a data.frame with information on all reconstructed population-level genomes (defaults to a variable named cI)
#               - "graph" - the network as graphNEL R-object (defaults to a variable named koGraph)
#               - "labCol" - a vector with colors to illustrate fold changes (defaults to a variable named hmcol)

#written by Anna Heintz-Buschart (July 2015)

library(BioNet)
library(RColorBrewer)
library(scales)

hmcol <- colorRampPalette(brewer.pal(9,"RdYlBu"))(512)[512:1]
compCol <- colorRampPalette(brewer.pal(11,"Spectral"))(109)[109:1]
cI <- readRDS("allClusterInfo.RDS")
koGraph <- readRDS("koGraph.RDS")

annaMustPrepModule <- function(nodes,diffRes,graph=koGraph,labCol=hmcol){
  modNW <- igraph.from.graphNEL(subNetwork(nodes,graph),name=T,weight=T)
  modLab <- merge(nodes,diffRes,by.x=1,by.y=0,all.x=T)$log2FoldChange 
  modLabC1 <- labCol[cut(c(modLab,-max(abs(modLab)),max(abs(modLab))),length(labCol))]
  modLabC2 <- labCol[length(labCol):1][cut(c(modLab,-max(abs(modLab)),max(abs(modLab))),length(labCol))]
  layo <- layout.fruchterman.reingold(modNW)
  return(list(modNW,modLabC1,modLabC2,layo))
}

annaMustPlotIndiModule <- function(nodes,nw,labCol,layout,indiv,outname,ome="T",silent=T,compc=compCol,clusterInf=cI){
  if(!all(ome %in% c("G","T","P"))){
    stop("invalid omic level, should be G,T or P")
  }else{
    if("G" %in% ome){
      indivExpr <- paste("perSampleNodes2Pop/",indiv,".DNAclusterPerNode.tsv",sep="")
      if(!file.exists(indivExpr)){
        stop("invalid individual, should be MX.X-VX")
      }else{
        nwPropG <- read.delim(indivExpr,stringsAsFactors=F)
        nwPropGmod <- merge(nodes,nwPropG,by.x=1,by.y=0,all.x=T)
        rownames(nwPropGmod) <- nwPropGmod[,1]
        nwPropGmod <- nwPropGmod[,-1]
        nwPropGmod[is.na(nwPropGmod)] <- 0
        nwPropGmod$absent <- ifelse(rowSums(nwPropGmod)==0,1,0)
        modCol <- c(sapply(colnames(nwPropGmod)[-ncol(nwPropGmod)],
                           function(x) compc[as.numeric(clusterInf$uniqueEss[clusterInf$cluster==x&clusterInf$sample==indiv])+1]),"grey")
        modAlph <- vector("numeric",length(modCol))
        first <- T
        for(i in 1:length(modAlph)){
          if(first) {
            first <- F
            modAlph[i] <- 1
          } else {
            if(cC==modCol[i]){
              modAlph[i] <- modAlph[i-1]*0.7
            }else{
              modAlph[i] <- 1
            }
          }
          cC <- modCol[i]
        }
        modCol <- list(alpha(modCol,modAlph))
        
        pdf(paste(outname,indiv,"metaG.pdf",sep="_"),width=7,height=7,pointsize=8)
        plot(nw,vertex.shape="pie",vertex.pie=as.list(as.data.frame(t(nwPropGmod))),
             vertex.size=4*log10(0.01+rowSums(nwPropGmod[,-1])), vertex.label.dist=.5,vertex.label.cex=.6,
             vertex.label.family="sans",vertex.pie.color=modCol,vertex.label.color=labCol,vertex.label.font=2,layout=layout)
        dev.off()
        if(!silent) plot(nw,vertex.shape="pie",vertex.pie=as.list(as.data.frame(t(nwPropGmod))),
                         vertex.size=4*log10(0.01+rowSums(nwPropGmod[,-1])), vertex.label.dist=.5,vertex.label.cex=.6,
                         vertex.label.family="sans",vertex.pie.color=modCol,vertex.label.color=labCol,vertex.label.font=2,
                         layout=layout)
      }
    }
    if("T" %in% ome){
      indivExpr <- paste("perSampleNodes2Pop/",indiv,".RNAclusterPerNode.tsv",sep="")
      if(!file.exists(indivExpr)){
        stop("invalid individual, should be MX.X-VX")
      }else{
        nwPropT <- read.delim(indivExpr,stringsAsFactors=F)
        nwPropTmod <- merge(nodes,nwPropT,by.x=1,by.y=0,all.x=T)
        rownames(nwPropTmod) <- nwPropTmod[,1]
        nwPropTmod <- nwPropTmod[,-1]
        nwPropTmod[is.na(nwPropTmod)] <- 0
        nwPropTmod$absent <- ifelse(rowSums(nwPropTmod)==0,1,0)
        modCol <- c(sapply(colnames(nwPropTmod)[-ncol(nwPropTmod)],
                           function(x) compc[as.numeric(clusterInf$uniqueEss[clusterInf$cluster==x&clusterInf$sample==indiv])+1]),"grey")
        modAlph <- vector("numeric",length(modCol))
        first <- T
        for(i in 1:length(modAlph)){
          if(first) {
            first <- F
            modAlph[i] <- 1
          } else {
            if(cC==modCol[i]){
              modAlph[i] <- modAlph[i-1]*0.7
            }else{
              modAlph[i] <- 1
            }
          }
          cC <- modCol[i]
        }
        modCol <- list(alpha(modCol,modAlph))
        
        pdf(paste(outname,indiv,"metaT.pdf",sep="_"),width=7,height=7,pointsize=8)
        plot(nw,vertex.shape="pie",vertex.pie=as.list(as.data.frame(t(nwPropTmod))),
             vertex.size=4*log10(0.01+rowSums(nwPropTmod[,-1])), vertex.label.dist=.5,vertex.label.cex=0.6,
             vertex.label.family="sans",vertex.pie.color=modCol,vertex.label.color=labCol,vertex.label.font=2,layout=layout)
        dev.off()
        if(!silent) plot(nw,vertex.shape="pie",vertex.pie=as.list(as.data.frame(t(nwPropTmod))),
                         vertex.size=4*log10(0.01+rowSums(nwPropTmod[,-1])), vertex.label.dist=.5,vertex.label.cex=.6,
                         vertex.label.family="sans",vertex.pie.color=modCol,vertex.label.color=labCol,vertex.label.font=2,
                         layout=layout)
      }
    }
    if("P" %in% ome){
      indivExpr <- paste("perSampleNodes2Pop/",indiv,".PROTclusterPerNode.tsv",sep="")
      if(!file.exists(indivExpr)){
        stop("invalid individual, should be MX.X-VX")
      }else{
        nwPropP <- read.delim(indivExpr,stringsAsFactors=F)
        nwPropPmod <- merge(nodes,nwPropP,by.x=1,by.y=0,all.x=T)
        rownames(nwPropPmod) <- nwPropPmod[,1]
        nwPropPmod <- nwPropPmod[,-1]
        nwPropPmod[is.na(nwPropPmod)] <- 0
        nwPropPmod$absent <- ifelse(rowSums(nwPropPmod)==0,1,0)
        modCol <- c(sapply(colnames(nwPropPmod)[-ncol(nwPropPmod)],
                           function(x) compc[as.numeric(clusterInf$uniqueEss[clusterInf$cluster==x&clusterInf$sample==indiv])+1]),"grey")
        modAlph <- vector("numeric",length(modCol))
        first <- T
        for(i in 1:length(modAlph)){
          if(first) {
            first <- F
            modAlph[i] <- 1
          } else {
            if(cC==modCol[i]){
              modAlph[i] <- modAlph[i-1]*0.7
            }else{
              modAlph[i] <- 1
            }
          }
          cC <- modCol[i]
        }
        modCol <- list(alpha(modCol,modAlph))
        
        pdf(paste(outname,indiv,"metaP.pdf",sep="_"),width=7,height=7,pointsize=8)
        plot(nw,vertex.shape="pie",vertex.pie=as.list(as.data.frame(t(nwPropPmod))),
             vertex.size=4*log10(0.01+rowSums(nwPropPmod[,-1])), vertex.label.dist=.5,vertex.label.cex=.6,
             vertex.label.family="sans",vertex.pie.color=modCol,vertex.label.color=labCol,vertex.label.font=2,layout=layout)
        dev.off()
        if(!silent) plot(nw,vertex.shape="pie",vertex.pie=as.list(as.data.frame(t(nwPropPmod))),
                         vertex.size=4*log10(0.01+rowSums(nwPropPmod[,-1])), vertex.label.dist=.5,vertex.label.cex=.6,
                         vertex.label.family="sans",vertex.pie.color=modCol,vertex.label.color=labCol,vertex.label.font=2,
                         layout=layout)
      }
    }
  }
}

plotModSuper <- function(nodes,diffRes,indiv,outname,groups=1,ome="T",silent=T,compc=compCol,clusterInf=cI,graph=koGraph,labCol=hmcol){
  if(any(groups>2)){
    stop("only 2 groups possible")
  }else{
    val <- annaMustPrepModule(nodes,diffRes,graph=koGraph,labCol=hmcol)
    if(length(groups)>1){
      for(i in unique(groups)){
        indi <- indiv[groups==i]
        for(ind in indi){
          annaMustPlotIndiModule(nodes,val[[1]],val[[1+i]],val[[4]],ind,outname,ome)
        }
      }
    }else{
      for(ind in indiv){
        annaMustPlotIndiModule(nodes,val[[1]],val[[1+groups]],val[[4]],ind,outname,ome)
      }
    }
  }
  return(list(val[[1]],val[[4]]))
}
