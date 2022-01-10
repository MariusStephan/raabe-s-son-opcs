#############   DE Wilcoxon testing for marker expression figures #############
# Version:      1
# Authour:      Marius Stephan
# Created:      27.09.2021
# Last edited:  27.09.2021

# Clear Workspace
rm(list=ls())

# Load Packages needed
library(Seurat)
library(tidyverse)
library(cowplot)
library(ggsignif)
library(ggrepel)

# Define local functions
graph.helper.sigAnnot <- function (p){
  if (is.na(p)){
    return(NA)
  } else if (p < 0.001) {
    return("***")
  } else if (p < 0.01) {
    return("**")
  } else if (p < 0.05) {
    return("*")
  } else if (p < 0.1) {
    return("p<0.1")
  } else {
    return("n.s.")
  }
}

# Load data
dm <- readRDS("OPCmerge_v09_OPCsubclusters.rds")

# Extract relevant features from normalised expression matrix
d <- as.matrix(dm@assays$SCT@data)[c("HES1","NES","PAX6","SOX3","MSI1","ASCL1","SOX9","ST8SIA1","ID4","CNP","PLP1","GALC","MBP","MOG","MOBP"),]
d <- data.frame(t(d))
d$cell <- row.names(d)
m <- dm@meta.data
m$cell <- row.names(m)

d <- left_join(d,select(m,cell,group,clusters),by="cell")
d <- gather(d, key = "gene", value = "expression",-cell,-group,-clusters)

## Compute Wilcoxon Rank-sum Tests for between sample comparisons with FDR adjustment
stats3 <- data.frame()
for (feat in unique(d$gene)){
  for (clust in unique(d$clusters)){
    levels <- unique(filter(d,gene == feat, clusters == clust)$group)
    
    if(length(levels)>1){
      w <- wilcox.test(filter(d, group == levels[1] & clusters == clust & gene == feat)$expression,filter(d, group == levels[2] & clusters == clust & gene == feat)$expression)
      w <- data.frame(W=unname(w[[1]]), p_val = w[[3]])
      w$n1 <- length(filter(d,gene == feat & clusters == clust & group == levels[1])$group)
      w$n2 <- length(filter(d,gene == feat & clusters == clust & group == levels[2])$group)
      w$Group1 <- levels[1]
      w$Group2 <- levels[2]
      
      if(length(levels)>2){
        w2 <- wilcox.test(filter(d,gene == feat & clusters == clust &group == levels[1])$expression,filter(d,gene == feat & clusters == clust &group == levels[3])$expression)
        w2 <- data.frame(W=unname(w2[[1]]), p_val = w2[[3]])
        w2$n1 <- length(filter(d,gene == feat & clusters == clust &group == levels[1])$group)
        w2$n2 <- length(filter(d,gene == feat & clusters == clust &group == levels[3])$group)
        w2$Group1 <- levels[1]
        w2$Group2 <- levels[3]
        
        w3 <- wilcox.test(filter(d,gene == feat, clusters == clust, group == levels[3])$expression,filter(d,gene == feat & clusters == clust &group == levels[2])$expression)
        w3 <- data.frame(W=unname(w3[[1]]), p_val = w3[[3]])
        w3$n1 <- length(filter(d,gene == feat & clusters == clust &group == levels[3])$group)
        w3$n2 <- length(filter(d,gene == feat & clusters == clust &group == levels[2])$group)
        w3$Group1 <- levels[3]
        w3$Group2 <- levels[2]
        
        w <- rbind(w,w2,w3)
      }
      w$Gene <- feat
      w$Cluster <- clust
      w$p_val_adj <- p.adjust(w$p_val)
      stats3 <- rbind(stats3,w)
    }
  }
}

stats3$p_val_adj_annot <- sapply(stats3$p_val_adj,graph.helper.sigAnnot)
write.table(stats3, "wilcox_sampleComparisons_Expression.tsv",quote = F,row.names = F,dec = ",",col.names = T,sep = "\t")


## Compute Wilcoxon Rank-sum Tests for between cluster comparisons with FDR adjustment
stats3 <- data.frame()
for (feat in unique(d$gene)){
  w <- wilcox.test(filter(d, gene == feat & clusters == "NPC")$expression, filter(d, gene == feat & clusters == "OPC")$expression)
  w <- data.frame(W=unname(w[[1]]), p_val = w[[3]])
  w$n1 <- length(filter(d, gene == feat & clusters == "NPC")$expression)
  w$n2 <- length(filter(d, gene == feat & clusters == "OPC")$expression)
  w$Group1 <- "NPC"
  w$Group2 <- "OPC"
  
  t <- wilcox.test(filter(d, gene == feat & clusters == "NPC")$expression, filter(d, gene == feat & clusters == "OPC")$expression)
  t <- data.frame(W=unname(t[[1]]), p_val = t[[3]])
  t$n1 <- length(filter(d, gene == feat & clusters == "NPC")$expression)
  t$n2 <- length(filter(d, gene == feat & clusters == "OPC")$expression)
  t$Group1 <- "NPC"
  t$Group2 <- "OPC"
  w <- rbind(w,t)
  
  t <- wilcox.test(filter(d, gene == feat & clusters == "OPC")$expression, filter(d, gene == feat & clusters == "OL1")$expression)
  t <- data.frame(W=unname(t[[1]]), p_val = t[[3]])
  t$n1 <- length(filter(d, gene == feat & clusters == "OPC")$expression)
  t$n2 <- length(filter(d, gene == feat & clusters == "OL1")$expression)
  t$Group1 <- "OPC"
  t$Group2 <- "OL1"
  w <- rbind(w,t)
  
  t <- wilcox.test(filter(d, gene == feat & clusters == "OL1")$expression, filter(d, gene == feat & clusters == "OL2")$expression)
  t <- data.frame(W=unname(t[[1]]), p_val = t[[3]])
  t$n1 <- length(filter(d, gene == feat & clusters == "OL1")$expression)
  t$n2 <- length(filter(d, gene == feat & clusters == "OL2")$expression)
  t$Group1 <- "OL1"
  t$Group2 <- "OL2"
  w <- rbind(w,t)
  
  w$Gene <- feat
  w$p_val_adj <- p.adjust(w$p_val)
  stats3 <- rbind(stats3,w)
}


stats3$p_val_adj_annot <- sapply(stats3$p_val_adj,graph.helper.sigAnnot)
write.table(stats3, "wilcox_clusterComparisons_Expression.tsv",quote = F,row.names = F,dec = ",",col.names = T,sep = "\t")
