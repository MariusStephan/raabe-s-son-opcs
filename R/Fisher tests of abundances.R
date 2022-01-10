#############   Chi-Squared Tests on Relevant Expressing Cell Abundances #############
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

# Local functions
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

### Import data
props <- read.csv("ClusterCell_Proportions_newSeurat.tsv",header = T,stringsAsFactors = F,sep="\t")

### Chi Squared Test to test for influence of sample on cluster abundances
chisq.test(select(spread(select(props,group,clusters,cluster_cells), key = group, value = cluster_cells, fill = 0),-clusters))

### Fisher's Exact Tests to test for influence of sample on cluster abundances
clusters <- unique(props$clusters)
props_fisher <- data.frame()
for (clust in clusters){
  
  d <- select(filter(props,clusters == clust),-percentage)
  
  d$rest <- d$total_cells - d$cluster_cells # Number of cells not in the respective cluster
  d <- rbind(select(d,-rest,-total_cells),data.frame(group=d$group,clusters="Rest",cluster_cells=d$rest))
  
  levels <- unique(d$group)
  
  if(length(levels)>1){
    
    t <- fisher.test(select(spread(filter(d, group %in% levels[c(1:2)]), key = group, value = cluster_cells, fill = 0),-clusters))
    f <- data.frame(odds_ratio = t$estimate, conf_int_lower = t$conf.int[1], conf_int_upper = t$conf.int[2],p_val = t$p.value, Group1 = levels[1], Group2 = levels[2])
    
    if(length(levels)>2){
      t <- fisher.test(select(spread(filter(d, group %in% levels[c(1,3)]), key = group, value = cluster_cells, fill = 0),-clusters))
      t <- data.frame(odds_ratio = t$estimate, conf_int_lower = t$conf.int[1], conf_int_upper = t$conf.int[2],p_val = t$p.value, Group1 = levels[1], Group2 = levels[3])
      f <- rbind(f,t)
      
      t <- fisher.test(select(spread(filter(d, group %in% levels[c(2,3)]), key = group, value = cluster_cells, fill = 0),-clusters))
      t <- data.frame(odds_ratio = t$estimate, conf_int_lower = t$conf.int[1], conf_int_upper = t$conf.int[2],p_val = t$p.value, Group1 = levels[2], Group2 = levels[3])
      f <- rbind(f,t)
    }
    f$cluster <- clust
    f$p_val_adj <- p.adjust(f$p_val)
    props_fisher <- rbind(props_fisher,f)
  }
}

### Fisher's Exact Test to test for influence of sample on the abundance of expressing cells in the respective cluster

### Load data
dm <- readRDS("OPCmerge_v09_OPCsubclusters.rds")

### Get expression level of genes of interest
features <- c("HES1","NES","PAX6","SOX3","MSI1","ASCL1","SOX9","ST8SIA1","ID4","CNP","PLP1","GALC","MBP","MOG","MOBP")
abundances <- as.data.frame(t(as.matrix(dm@assays$RNA@counts[which(row.names(dm@assays$RNA@counts) %in% features),])))

### Compute abundance from expression levels
abundances$cell <- row.names(abundances)
abundances <- gather(abundances,key = gene,value = expression,-cell)
temp <- dm@meta.data
temp$cell <- row.names(temp)
temp <- select(temp,cell,group,clusters)
abundances <- left_join(abundances,temp,by = "cell")
abundances %>% group_by(gene,group,clusters) %>% summarise(expressing = sum(expression > 0), nonExpressing = sum(expression == 0)) -> abundances
abundances <- ungroup(abundances)


### Conduct test for every gene of interest and cluster (between samples)
abundances_fisher <- data.frame()
for (feat in features){
  for (clust in clusters){
    
    d <- filter(abundances, gene == feat, clusters == clust)
    
    levels <- unique(d$group)
    
    if(length(levels)>1){
      
      t <- fisher.test(select(filter(d,group %in% levels[c(1:2)]),expressing,nonExpressing))
      f <- data.frame(odds_ratio = t$estimate, conf_int_lower = t$conf.int[1], conf_int_upper = t$conf.int[2],p_val = t$p.value, Group1 = levels[1], Group2 = levels[2])
      
      if(length(levels)>2){
        t <- fisher.test(select(filter(d,group %in% levels[c(1,3)]),expressing,nonExpressing))
        t <- data.frame(odds_ratio = t$estimate, conf_int_lower = t$conf.int[1], conf_int_upper = t$conf.int[2],p_val = t$p.value, Group1 = levels[1], Group2 = levels[3])
        f <- rbind(f,t)
        
        t <- fisher.test(select(filter(d,group %in% levels[c(2:3)]),expressing,nonExpressing))
        t <- data.frame(odds_ratio = t$estimate, conf_int_lower = t$conf.int[1], conf_int_upper = t$conf.int[2],p_val = t$p.value, Group1 = levels[2], Group2 = levels[3])
        f <- rbind(f,t)
      }
      f$gene <- feat
      f$cluster <- clust
      f$p_val_adj <- p.adjust(f$p_val)
      abundances_fisher <- rbind(abundances_fisher,f)
    }
  }
}

### family-wise comparisons between clusters
fisherExpression_family <- data.frame()
for (feat in features){
  d <- filter(abundances, gene == feat)
  d %>% group_by(group) %>% summarise(expressing = sum(expressing),nonExpressing = sum(nonExpressing)) -> d
  t <- fisher.test(select(d,expressing,nonExpressing))
  fisherExpression_family <- rbind(fisherExpression_family, data.frame(p_val = t$p.value, gene = feat))
}

### pair-wise comparisons between clusters (with FDR adjustment)
fisherExpression_pair <- data.frame()
for (feat in features){
  d <- filter(abundances, gene == feat)
  d %>% group_by(group) %>% summarise(expressing = sum(expressing),nonExpressing = sum(nonExpressing)) -> d
  levels <- unique(d$group)
  
  if(length(levels)>1){
    
    t <- fisher.test(select(filter(d,group %in% levels[c(1:2)]),expressing,nonExpressing))
    f <- data.frame(odds_ratio = t$estimate, conf_int_lower = t$conf.int[1], conf_int_upper = t$conf.int[2],p_val = t$p.value, Group1 = levels[1], Group2 = levels[2])
    
    if(length(levels)>2){
      t <- fisher.test(select(filter(d,group %in% levels[c(1,3)]),expressing,nonExpressing))
      t <- data.frame(odds_ratio = t$estimate, conf_int_lower = t$conf.int[1], conf_int_upper = t$conf.int[2],p_val = t$p.value, Group1 = levels[1], Group2 = levels[3])
      f <- rbind(f,t)
      
      t <- fisher.test(select(filter(d,group %in% levels[c(2:3)]),expressing,nonExpressing))
      t <- data.frame(odds_ratio = t$estimate, conf_int_lower = t$conf.int[1], conf_int_upper = t$conf.int[2],p_val = t$p.value, Group1 = levels[2], Group2 = levels[3])
      f <- rbind(f,t)
    }
    f$gene <- feat
    f$p_val_adj <- p.adjust(f$p_val)
    fisherExpression_pair <- rbind(fisherExpression_pair,f)
    
  }
}

### Add annotation of statistical significance
props_fisher$p_val_adj_annot <- sapply(props_fisher$p_val_adj,graph.helper.sigAnnot)
abundances_fisher$p_val_adj_annot <- sapply(abundances_fisher$p_val_adj,graph.helper.sigAnnot)

fisherExpression_family$p_val_annot <- sapply(fisherExpression_family$p_val,graph.helper.sigAnnot)
fisherExpression_pair$p_val_adj_annot <- sapply(fisherExpression_pair$p_val_adj,graph.helper.sigAnnot)

### Export results
write.table(props_fisher, "Fisher_clusterDistribution.tsv",quote = F,row.names = F,dec = ",",col.names = T,sep = "\t")
write.table(abundances_fisher, "Fisher_expressionAbundance.tsv",quote = F,row.names = F,dec = ",",col.names = T,sep = "\t")

write.table(fisherExpression_family, "Fisher_expressionAbundance_familywise.tsv",quote = F,row.names = F,dec = ",",col.names = T,sep = "\t")
write.table(fisherExpression_pair, "Fisher_expressionAbundance_pairwise.tsv",quote = F,row.names = F,dec = ",",col.names = T,sep = "\t")
