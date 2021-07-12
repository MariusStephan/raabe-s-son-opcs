#############   Export Seurat dataset for import in Velocyto (Python)  #############
# Version:      03 (OPC subclusters)
# Authour:      Marius Stephan
# Created:      29.03.2021
# Last edited:  18.05.2021

# Clear Workspace
rm(list=ls())

# Load Packages needed
library(Seurat)
library(scales)
library(tidyverse)
library(SeuratDisk)
library(scales)

# Load dataset from R dump
dm <- readRDS("C:/Users/Marius/MNB QSync/Exchange/Marius/RStudio Projects/2021/2021-01-19 Florian NPC/OPCmerge_v09_OPCsubclusters.rds")
met <- as.numeric(dm@meta.data$seurat_clusters)
names(met) <- names(dm@active.ident)

### Extract Cell IDs, Cluster IDs and UMAP individually

# Extract Cell IDs and save them to file
write.csv(Cells(dm), file = "OPCmerged_cellIDs.csv", row.names = FALSE)

# Extract UMAP coordinates and save them to file
write.csv(Embeddings(dm, reduction = "umap"), file = "OPCmerged_UMAP.csv")

# Extract Cluster IDs and save them to file
write.csv(dm@active.ident, file = "OPCmerged_clusters_OPCsplitClusters.csv")
write.csv(names(dm@active.ident)[dm@active.ident != 'Neuronal'], file = "OPCmerged_cellIDs_noNeuronalCells.csv")
write.csv(names(dm@active.ident)[dm@active.ident %in% c('Neuronal','NPC')], file = "OPCmerged_cellIDs_NeuronalLineage.csv")

# Extract Cluster numbers and save them to file
write.csv(met, file = "OPCmerged_clusterNumbers.csv")

# Get cluster colors
identities <- levels(dm@meta.data$clusters)
color_palette <- hue_pal()(length(identities))
names(color_palette) <- identities
write.csv(color_palette, file = "OPCmerged_clusterPalette.csv")

# Get updated cluster colors (with OPC1 OPC2 division)
identities <- levels(dm@meta.data$clusters2)
color_palette <- c("#fcc8c4","#944641","#A3A500","#00BF7D","#00B0F6","#E76BF3")
names(color_palette) <- identities
write.csv(color_palette, file = "OPCmerged_clusterPalette_splitClusters.csv")

color_palette <- data.frame(clusters =identities, colors = color_palette, stringsAsFactors = F)

color_palette <- left_join(select(dm@meta.data, clusters),color_palette)
row.names(color_palette) <- names(dm@active.ident)
color_palette <- select(color_palette, -clusters)
write.csv(color_palette, file = "OPCmerged_clusterColors.csv")

# Get samples colors
rawColors <- as.character(as.hexmode(col2rgb(c("gray","light blue","dark blue"))))
groupColors <- c()
for(i in 1:ncol(rawColors)){
  groupColors[i] <- paste0("#",rawColors[1,i],rawColors[2,i],rawColors[3,i])
}
groupColors <- data.frame(x=c("Pre","S","SON"),colors=groupColors)
groupPalette <-  toupper(groupColors$colors)
names(groupPalette) <- groupColors$x
write.csv(groupPalette, file = "OPCmerged_groupPalette.csv")

# Get sample IDs and add Cell ID
d <- dm@meta.data$group
names(d) <- names(dm@active.ident)
write.csv(d, file = "OPCmerged_groups.csv")
