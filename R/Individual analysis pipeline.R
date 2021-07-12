#############   NPC/OPC Individual Sample Analysis  #############
# Version:      01
# Authour:      Marius Stephan
# Created:      06.05.2021
# Last edited:  06.05.2021

# Clear Workspace
rm(list=ls())

# Load Packages needed
library(Seurat)
library(tidyverse)
library(ggrepel)
library(cowplot)

### Local functions
VolcanoPlot <- function (data,log2FC_threshold = 0.25, p_val_threshold = 0.05, adjust_p_val = FALSE, max_labeled_hits = 20, x_symmetry = TRUE){
  if("avg_log2FC" %in% colnames(data) & "p_val" %in% colnames(data)){
    if (adjust_p_val){
      if ("p_val_adj" %in% colnames(data)){
        data$logp <- -log10(data$p_val_adj)
      } else {
        errorCondition("Dataset is not formatted as expected by this function!")      
      }
      
    } else {
      data$logp <- -log10(data$p_val)
    }
    
    pos_hitlist <- filter(data, logp >= -log10(p_val_threshold) & avg_log2FC > log2FC_threshold)
    neg_hitlist <- filter(data, logp >= -log10(p_val_threshold) & avg_log2FC < -log2FC_threshold)
    rest <- filter(data, logp < -log10(p_val_threshold) | abs(avg_log2FC) < log2FC_threshold)
    
    if(nrow(pos_hitlist) > 0){
      x.limits <-  max(pos_hitlist$avg_log2FC)
    } else {
      x.limits <-  max(rest$avg_log2FC)
    }
    
    if(nrow(neg_hitlist) > 0){
      x.limits <-  c(min(neg_hitlist$avg_log2FC),x.limits)
    } else{
      x.limits <-  c(min(rest$avg_log2FC), x.limits)
    }
    
    if(x_symmetry){
      x.limits <- c(-max(abs(x.limits)),max(abs(x.limits)))
    }
    
    if(nrow(pos_hitlist) > max_labeled_hits) {
      pos_hits2label <- slice(arrange(pos_hitlist,desc(avg_log2FC)),1:max_labeled_hits)
    } else {
      pos_hits2label <- pos_hitlist
    }
    
    if(nrow(neg_hitlist) > max_labeled_hits) {
      neg_hits2label <- slice(arrange(neg_hitlist,avg_log2FC),1:max_labeled_hits)
    } else {
      neg_hits2label <- neg_hitlist
    }
    
    plot <- ggplot(mapping = aes(x = avg_log2FC, y = logp))+
      geom_point(data=rest, color = "gray",size = 3)+
      geom_point(data=neg_hitlist, color = "dark blue",size = 3)+
      geom_point(data=pos_hitlist, color = "dark red",size = 3)+
      geom_hline(yintercept = -log10(p_val_threshold),lty=2, colour = "gray",size = 1.2)+
      geom_vline(xintercept = log2FC_threshold,lty=2, colour = "gray",size = 1.2)+
      geom_vline(xintercept = -log2FC_threshold,lty=2, colour = "gray",size = 1.2)+
      geom_text_repel(data = neg_hits2label, aes(label = gene), colour = " dark blue",force = 5,point.padding = 1)+
      geom_text_repel(data = pos_hits2label, aes(label = gene), colour = "dark red",force = 5, point.padding = 1)+
      theme_cowplot()+
      xlim(x.limits)+
      xlab('log2FC(A/B)')+
      ylab('-log10(p)')
    
    return(plot)
    
  } else {
    errorCondition("Dataset is not formatted as expected by this function!")
  }
}

### Load R Data Dump to get cell IDs
dm <- readRDS("C:/Users/Marius/MNB QSync/Exchange/Marius/RStudio Projects/2021/2021-01-19 Florian NPC/OPCmerge_v09_OPCsubclusters.rds")

### Extract cell list, UMAP coordinates and clusters
cells <- names(dm@active.ident)
rm(dm)
### Read count matrices

# Import DGE data from Drop-Seq pipeline (assuming 2000 cells went into DGE!!!)
d1 <- read.table(file = "Raw/NPC_S1_Construct_DGE_2000cells.txt", header = TRUE, row.names = 1, colClasses =c("character", rep("numeric", 2000)))
d2 <- read.table(file = "Raw/OPC1_S_S5_Construct_DGE_2000cells.txt", header = TRUE, row.names = 1, colClasses =c("character", rep("numeric", 2000)))
d3 <- read.table(file = "Raw/OPC2_SON_S6_Construct_DGE_2000cells.txt", header = TRUE, row.names = 1, colClasses =c("character", rep("numeric", 2000)))

## Filter count matrices for TFs
# We don't want the construct expression to influence downstream analyses, but look at the endogenous gene expression only.
d1 <- d1[row.names(d1) != "LTR",]
d2 <- d2[row.names(d2) != "LTR",]
d3 <- d3[row.names(d3) != "LTR",]


# Add sample ID to cell ID
colnames(d1) <- paste0("Pre_",colnames(d1))
colnames(d2) <- paste0("S_",colnames(d2))
colnames(d3) <- paste0("SON_",colnames(d3))

# Initialize the Seurat object with the raw (non-normalized data).
# Keep all genes expressed in >= 3 cells. Keep all cells with at least 200 detected genes
d1 <- CreateSeuratObject(counts = d1, min.cells = 3, project = "OPC")
d2 <- CreateSeuratObject(counts = d2, min.cells = 3, project = "OPC")
d3 <- CreateSeuratObject(counts = d3, min.cells = 3, project = "OPC")

# Add group identifier
d1@meta.data$group <- "NI Day+0"
d2@meta.data$group <- "S Day+15"
d3@meta.data$group <- "SON Day+15"

### Filter for cells used in merged dataset. We found that this step doesn't have a major influence on the results, but cells were already filtered in each sample individually in the merged pipeline. We've chosen to do it this way for consistency.
d1 <- subset(d1, cells = cells)
d2 <- subset(d2, cells = cells)
d3 <- subset(d3, cells = cells)

d1[['percent.mt']] <- PercentageFeatureSet(d1, pattern = "^MT-")
d2[['percent.mt']] <- PercentageFeatureSet(d2, pattern = "^MT-")
d3[['percent.mt']] <- PercentageFeatureSet(d3, pattern = "^MT-")

d1 <- SCTransform(d1, vars.to.regress = "percent.mt")
d2 <- SCTransform(d2, vars.to.regress = "percent.mt")
d3 <- SCTransform(d3, vars.to.regress = "percent.mt")

d1 <- FindVariableFeatures(d1, nfeatures = 3000)
d1 <- RunPCA(d1)

d2 <- FindVariableFeatures(d2, nfeatures = 3000)
d2 <- RunPCA(d2)

d3 <- FindVariableFeatures(d3, nfeatures = 3000)
d3 <- RunPCA(d3)

print(ElbowPlot(d1))
print(ElbowPlot(d2))
print(ElbowPlot(d3))

### I've found the best method to choose the right number of PCs to consider is to test the range at the joint of the elbow plot and watch how the result changes.
# You will often find that the differences are minor except for directionality (axes are often flipped or inverted).

# pdf("UMAP_Individual_Samples.pdf")
# for (pc in 10:20) {
# d1 <- FindNeighbors(d1, dims = 1:pc)
# d1 <- FindClusters(d1, resolution = 0.16)
# d1 <- RunUMAP(d1,dims = 1:pc)
# print(DimPlot(d1,reduction = "umap", pt.size = 2)+ggtitle(paste0("Pre ",pc," PCs")))
# }
# for (pc in 10:20) {
#   d2 <- FindNeighbors(d2, dims = 1:pc)
#   d2 <- FindClusters(d2, resolution = 0.16)
#   d2 <- RunUMAP(d2,dims = 1:pc)
#   print(DimPlot(d2,reduction = "umap", pt.size = 2)+ggtitle(paste0("S ",pc," PCs")))
# }
# for (pc in 10:20) {
#   d3 <- FindNeighbors(d3, dims = 1:pc)
#   d3 <- FindClusters(d3, resolution = 0.16)
#   d3 <- RunUMAP(d3,dims = 1:pc)
#   print(DimPlot(d3,reduction = "umap", pt.size = 2)+ggtitle(paste0("SON ",pc," PCs")))
# }
# dev.off()

### Run kNN and UMAP
d1 <- FindNeighbors(d1, dims = 1:14)
d1 <- RunUMAP(d1,dims = 1:14)
d2 <- FindNeighbors(d2, dims = 1:20)
d2 <- RunUMAP(d2,dims = 1:20)
d3 <- FindNeighbors(d3, dims = 1:19)
d3 <- RunUMAP(d3,dims = 1:19)

### Test hierarchical clustering resolutions.
# I increase the the resolution until all properly discernible clusters visible in the UMAP embedding are separated and stop at that point to avoid overclustering.

# pdf("kNN_Clustering_Individual_Samples.pdf")
# for (res in seq(0.05,0.5,0.05)) {
#   dm <- FindClusters(d1, resolution = res)
#   print(DimPlot(dm,reduction = "umap", pt.size = 2)+ggtitle(paste0("Pre ",res," Resolution")))
# }
# for (res in seq(0.05,0.5,0.05)) {
#   dm <- FindClusters(d2, resolution = res)
#   print(DimPlot(dm,reduction = "umap", pt.size = 2)+ggtitle(paste0("S ",res," Resolution")))
# }
# for (res in seq(0.05,0.5,0.05)) {
#   dm <- FindClusters(d3, resolution = res)
#   print(DimPlot(dm,reduction = "umap", pt.size = 2)+ggtitle(paste0("SON ",res," Resolution")))
# }
# dev.off()

### Run clustering
d1 <- FindClusters(d1, resolution = 0.3)
d2 <- FindClusters(d2, resolution = 0.3)
d3 <- FindClusters(d3, resolution = 0.35)

### Annotate clusters
new.cluster.ids <- c("NPC","Neuronal","Unspecified")
names(new.cluster.ids) <- levels(d1)
d1 <- RenameIdents(d1, new.cluster.ids)
levels(d1) <- c("Neuronal","NPC","Unspecified")
d1@meta.data$clusters<- unname(d1@active.ident)

new.cluster.ids <- c("OPC1","OPC2","Neuronal","OL")
names(new.cluster.ids) <- levels(d2)
d2 <- RenameIdents(d2, new.cluster.ids)
levels(d2) <- c("Neuronal","OPC1","OPC2","OL")
d2@meta.data$clusters<- unname(d2@active.ident)

new.cluster.ids <- c("OPC","OL1","OL2","Neuronal")
names(new.cluster.ids) <- levels(d3)
d3 <- RenameIdents(d3, new.cluster.ids)
levels(d3) <- c("Neuronal","OPC","OL1","OL2")
d3@meta.data$clusters<- unname(d3@active.ident)

### Create UMAP plots
pdf("UMAP_DimPlots_Individual_Samples.pdf",height=6,width=8)
print(DimPlot(d1,reduction = "umap", pt.size = 2,label = T))
print(DimPlot(d2,reduction = "umap", pt.size = 2,label = T))
print(DimPlot(d3,reduction = "umap", pt.size = 2,label = T))
dev.off()

### Look for DEGs
markers <- FindMarkers(d2,ident.1 = "OPC2", ident.2 = "OPC1", assay = "SCT",logfc.threshold = 0)
markers$Cluster1 <- "OPC2"
markers$Cluster2 <- "OPC1"
markers$gene <- row.names(markers)
#markers <- filter(markers, gene %in% tf)
markers <- arrange(markers,desc(avg_log2FC),p_val_adj)
write.table(markers,"SCT_DGE_individualAnalysis_S_OPC2vsOPC1_unfiltered_forGSEA.tsv",quote = F,row.names = F,col.names = T,sep = "\t")

VolcanoPlot(markers,log2FC_threshold = 1)+ggtitle("S Individual Analysis OPC2 vs. OPC1")

markers <- FindMarkers(d3,ident.1 = "OL2", ident.2 = "OL1", assay = "SCT",logfc.threshold = 0)
markers$Cluster1 <- "OL2"
markers$Cluster2 <- "OL1"
markers$gene <- row.names(markers)
#markers <- filter(markers, gene %in% tf)
markers <- arrange(markers,desc(avg_log2FC),p_val_adj)
write.table(markers,"SCT_DGE_individualAnalysis_SON_OL2vsOL1_unfiltered_forGSEA.tsv",quote = F,row.names = F,col.names = T,sep = "\t")

pdf("VolcanoPlot_SCT_DGE_individualAnalysis_SON_OL2vsOL1.pdf",width = 8,height = 6)
print(VolcanoPlot(markers,log2FC_threshold = 1)+ggtitle("S Individual Analysis OPC2 vs. OPC1"))
dev.off()
