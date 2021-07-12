#############   SureCell WTA 3' scRNA-Seq Integrated Data Analysis  #############
# Version:      04 (S SON Construct added)
# Authour:      Marius Stephan
# Created:      19.01.2021
# Last edited:  06.07.2021

# Clear Workspace
rm(list=ls())

# Load Packages needed
library(Seurat)
library(tidyverse)
library(patchwork)
library(Rmagic)
library(cowplot)
library(ggsignif)
library(pheatmap)
library(ggrepel)

# Load user-defined functions
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

#####################################################
#  saveRDS(object = dm,file = "C:/Users/Marius/MNB QSync/Exchange/Marius/RStudio Projects/2021/2021-01-19 Florian NPC/OPCmerge_v09_OPCsubclusters.rds")
# dm <- readRDS("C:/Users/Marius/MNB QSync/Exchange/Marius/RStudio Projects/2021/2021-01-19 Florian NPC/OPCmerge_v09_OPCsubclusters.rds")
# dm <- readRDS("C:/Users/Marius/MNB QSync/Exchange/Marius/RStudio Projects/2021/2021-01-19 Florian NPC/OPCmerge_v08_latentTime.rds")
#####################################################


### Data import and preprocessing
############

# Import DGE data from Drop-Seq pipeline (assuming 2000 cells went into DGE!!!)
 d1 <- read.table(file = "Raw/NPC_S1_Construct_DGE_2000cells.txt.gz", header = TRUE, row.names = 1, colClasses =c("character", rep("numeric", 2000)))
 d2 <- read.table(file = "Raw/OPC1_S_S5_Construct_DGE_2000cells.txt.gz", header = TRUE, row.names = 1, colClasses =c("character", rep("numeric", 2000)))
 d3 <- read.table(file = "Raw/OPC2_SON_S6_Construct_DGE_2000cells.txt.gz", header = TRUE, row.names = 1, colClasses =c("character", rep("numeric", 2000)))

# Add sample ID to cell ID
 colnames(d1) <- paste0("Pre_",colnames(d1))
 colnames(d2) <- paste0("S_",colnames(d2))
 colnames(d3) <- paste0("SON_",colnames(d3))

# Filter count matrices for TFs
# We want to avoid the construct reads to influence downstream gene expression analyses. Note that the dataset in the paper was acquired without mapping of construct reads at first, latter being added downstream in the process.
# This does not have major influence on the results of this study. However, even slightest changes in the count matrix such as column switches may cause the axes of the PCA and UMAP results to flip.
d1 <- d1[row.names(d1) != "LTR",]
d2 <- d2[row.names(d2) != "LTR",]
d3 <- d3[row.names(d3) != "LTR",]
 
# Create the Seurat object
# Keep all genes expressed in >= 3 cells (~0.1% of the data). Keep all cells with at least 200 detected genes
 d1 <- CreateSeuratObject(counts = d1, min.cells = 3, min.features = 200, project = "OPC")
 d2 <- CreateSeuratObject(counts = d2, min.cells = 3, min.features = 200, project = "OPC")
 d3 <- CreateSeuratObject(counts = d3, min.cells = 3, min.features = 200, project = "OPC")

# Add group identifier
 d1@meta.data$group <- "Pre"
 d2@meta.data$group <- "S"
 d3@meta.data$group <- "SON"

### QC, merge and normalization
#############
 
# Find abundance of mitochondrial transcipt (High abundance marks apoptotic/damaged cells)
 d1[['percent.mt']] <- PercentageFeatureSet(d1, pattern = "^MT-")
 d2[['percent.mt']] <- PercentageFeatureSet(d2, pattern = "^MT-")
 d3[['percent.mt']] <- PercentageFeatureSet(d3, pattern = "^MT-")

# QC Violin Plots
print(VlnPlot(d1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3))
print(VlnPlot(d2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3))
print(VlnPlot(d3, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3))

#### Filter Data with 'recommended' settings
d1 <- subset(d1, subset = nCount_RNA > 1700 & nFeature_RNA > 1000 & nFeature_RNA < 6000 & percent.mt < 5) 
print(paste0(length(d1$nFeature_RNA)," cells out of total input fulfill the quality criteria."))

d2 <- subset(d2, subset = nCount_RNA > 1000 & nFeature_RNA > 1000 & nFeature_RNA < 6000 & percent.mt < 5) 
print(paste0(length(d2$nFeature_RNA)," cells out of total input fulfill the quality criteria."))

d3 <- subset(d3, subset = nCount_RNA > 1700 & nFeature_RNA > 1000 & nFeature_RNA < 6000 & percent.mt < 5)
print(paste0(length(d3$nFeature_RNA)," cells out of total input fulfill the quality criteria."))


### Put datasets together in a list
dList <- list(d1,d2,d3)
rm(d1,d2,d3)

### Merge instead of integration
dm <- merge(x = dList[[1]],y = dList[2:3], merge.data = F)
 
#### Normalize and scale data with sctransform
dm <- SCTransform(dm,vars.to.regress =  c("percent.mt")) 


### Dimension reduction, clustering and UMAP embedding
#############

#### First Dimension Reduction Step (PCA)
 dm <- FindVariableFeatures(dm, nfeatures = 3000)
 dm <- RunPCA(dm)

# Check PCA 
print(DimPlot(dm, reduction = "pca",dims = c(1,2)))
print(DimPlot(dm, reduction = "pca",dims = c(2,3)))

# Elbow plot
print(ElbowPlot(dm))

### Clustering
# Note the resolution value is determined in part by trial and error.
# One might want to revisit this section after having the clusters visualized in a UMAP DimPlot

### Test PC number in the range of the knee plot "joint" 
 pdf("UMAP.pdf")
 for (pc in 10:20) {
 dm <- FindNeighbors(dm, dims = 1:pc)
 dm <- FindClusters(dm, resolution = 0.16)
 dm <- RunUMAP(dm,dims = 1:pc)
 print(DimPlot(dm,reduction = "umap",group.by = "group", pt.size = 2)+ggtitle(paste0(pc," PCs")))
 print(DimPlot(dm,reduction = "umap",split.by = "group", pt.size = 2)+ggtitle(paste0(pc," PCs")))
 }
 dev.off()

 ### Run SNN clustering and UMAP embedding with appropriate PC number
 pc <- 15
 dm <- FindNeighbors(dm, dims = 1:pc)
 dm <- FindClusters(dm, resolution = 0.23)
 dm <- RunUMAP(dm,dims = 1:pc)

 ### Annotate clusters
 new.cluster.ids <- c("iOPC","NPC","Neuronal","iOL1","iOL2","unspecified") # unspecified cells show high expression of keratins and other epithelial markers. They might constitute an contamination or maybe a "cancerous" epithelial cell population
 names(new.cluster.ids) <- levels(dm)
 dm <- RenameIdents(dm, new.cluster.ids)
 dm$clusters <- Idents(dm)
 print(DimPlot(dm,reduction = "umap", pt.size = 2)+ggtitle(paste0(pc," PCs")))
 VlnPlot(dm,c("MBP"))+theme(legend.position = "right")
 
 ### Remove potential contamination
 dm <- subset(dm, idents = c("iOPC","NPC","Neuronal","iOL1","iOL2"))
 pc <- 15
 dm <- RunUMAP(dm,dims = 1:pc) # UMAP is rerun to readjust latent variables for the effect of the (now missing) contaminant cells
 new.cluster.ids <- factor(dm@meta.data$clusters,levels = c("Neuronal","NPC","iOPC","iOL1","iOL2"))
 names(new.cluster.ids) <- names(dm@active.ident)
 Idents(dm) <- new.cluster.ids
 VlnPlot(dm,c("MBP"))+theme(legend.position = "right")
 
 ### UMAP embedding dimplots
 print(DimPlot(dm,repel = T, reduction = "umap", group.by =  "group", pt.size = 2,shape.by = "clusters"))
 print(DimPlot(dm,reduction = "umap", pt.size = 2)+scale_color_manual(values=c("#00BF7D", "#A3A500", "#F8766D", "#E76BF3", "#00B0F6")))
 print(DimPlot(dm,repel = T, reduction = "umap", split.by =  "group", pt.size = 2)+scale_color_manual(values=c("#00BF7D", "#A3A500", "#F8766D", "#E76BF3", "#00B0F6")))
 print(DimPlot(dm,repel = T, reduction = "umap", group.by =  "group", pt.size = 2,shape.by = "clusters")+scale_color_manual(values = c("light gray","light blue","dark blue")))
 print(DimPlot(dm,repel = T, reduction = "umap", group.by =  "group", pt.size = 2,shape.by = "clusters")+scale_color_manual(values = c("light gray","light blue","dark blue")))
 print(DimPlot(dm,repel = T, reduction = "umap",label = TRUE, pt.size = 2,label.size = 5)+scale_color_manual(values=c("#00BF7D", "#A3A500", "#F8766D", "#00B0F6", "#E76BF3")))
 
 ### Add construct expression to dataset
 ###############
 
 ### Flag was not detected because of its position far from the downstream LTR, which serves as polyA site in the construct.
 # For this step to work, you need to skip the LTR deletion step. In the analysis for the paper, we've added these columns to the Seurat object previously saved in an R data dump (see below).
 construct_expression <- dm@assays$SCT@data["LTR",]
 construct_counts <- dm@assays$SCT@counts["LTR",]
 FeaturePlot(dm, repel = T, split.by = "group",sort.cell = T,slot = "counts" , features = "LTR",reduction = "umap",cols = c("light blue","dark red"),ncol = 1,pt.size = 1.5,label = T,label.size = 5,keep.scale = 'all')+theme(legend.position = "right")
 FeaturePlot(dm, repel = T, split.by = "group",sort.cell = T, features = "LTR",reduction = "umap",cols = c("light blue","dark red"),ncol = 1,pt.size = 1.5,label = T,label.size = 5,keep.scale = 'all')+theme(legend.position = "right")
 
 ### Add construct expression to meta data
 meta <- dm@meta.data
 meta$CB <- row.names(meta)
 construct_expression <- data.frame(CB=names(construct_counts),construct_counts=construct_counts,construct_normExpression=construct_expression,stringsAsFactors = F)
 meta <- left_join(meta,construct_expression)
 row.names(meta) <- meta$CB
 meta <- select(meta, -CB)
 dm@meta.data <- meta ### You might add this to previously saved RData dump (rds file)
 
pdf("ConstructExpression_UMAP.pdf",width = 6, height = 5)
 FeaturePlot(dm, repel = T,sort.cell = T, features = "construct_normExpression",reduction = "umap",cols = c("light blue","dark red"),ncol = 1,pt.size = 1.5,label = T,label.size = 5,keep.scale = 'all')+theme(legend.position = "right")
dev.off()

pdf("ConstructExpression_UMAP_splitBySample.pdf",width = 16, height = 5)
 FeaturePlot(dm, repel = T,sort.cell = T,split.by = "group", features = "construct_normExpression",reduction = "umap",cols = c("light blue","dark red"),ncol = 1,pt.size = 1.5,label = T,label.size = 5,keep.scale = 'all')+theme(legend.position = "right")
dev.off() 

pdf("ConstructExpression_VlnPlot.pdf",width = 8, height = 6)
VlnPlot(dm, pt.size = 0, features = "construct_normExpression",split.by = "group", cols = c("light gray", "light blue", "dark blue"))
VlnPlot(subset(dm, group != "NI Day+0"), pt.size = 0, features = "construct_normExpression",split.by = "group", cols = c("light blue", "dark blue"))
dev.off() 

### Add gene-shared latent time and velocity pseudotime
############

 ### Add similarity-based diffusion pseudo time / gene-shared latent time to meta data
 time <- read.delim("OPCmerge_pseudotime.tsv")
 time <- rename(time, CB = 'X')

 meta <- dm@meta.data
 meta$CB <- row.names(meta)
 meta <- left_join(meta,time)
 row.names(meta) <- meta$CB
 meta <- select(meta, -CB)
 dm@meta.data <- meta
 
 FeaturePlot(dm, repel = T,sort.cell = T, features = "velocity_pseudotime",reduction = "umap",cols = c("light blue","dark red"),ncol = 1,pt.size = 1.5,label = T,label.size = 5,keep.scale = 'all')+theme(legend.position = "right")
 FeaturePlot(dm, repel = T,sort.cell = T, features = "latent_time",reduction = "umap",cols = c("light blue","dark red"),ncol = 1,pt.size = 1.5,label = T,label.size = 5,keep.scale = 'all')+theme(legend.position = "right")
 
 pdf('ConstructExpression_Vlnplots_v02.pdf',height = 5,width = 6)
 print(VlnPlot(dm,"construct_normExpression",group.by = "group",cols = c("light gray","light blue","dark blue"),pt.size = 0)+theme(legend.position = "right"))
 print(VlnPlot(dm,"construct_normExpression",split.by = "group",cols = c("light gray","light blue","dark blue"),pt.size = 0)+theme(legend.position = "right"))
 print(VlnPlot(dm,"construct_normExpression",pt.size = 0)+theme(legend.position = "right")+scale_fill_manual(values=c("#00BF7D", "#A3A500", "#F8766D", "#E76BF3", "#00B0F6")))
 dev.off()

 ### Get cluster markers, abundances & average gene expression for later reference
 ############
 
### Cluster Proportions
 dm@meta.data %>% group_by(group) %>% summarise(total_cells = n()) -> props
 dm@meta.data %>% group_by(group, clusters) %>% summarise(cluster_cells = n()) -> props.cluster
 props <- left_join(props.cluster,props)
 props$percentage <- props$cluster_cells/props$total_cells*100
 write.table(props,"Project/ClusterCell_Proportions_newSeurat.tsv",sep="\t",quote = F,row.names = F,col.names = T)

 #### Find Cluster Markers
 dm.markers <- FindAllMarkers(dm, only.pos = T, min.pct = 0.1, logfc.threshold = 0.25,assay = "SCT")
 dm.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)
 write.table(dm.markers,"Cluster_PosMarkers_SCTransform_newSeurat.txt",row.names = F,col.names = T,quote = F,sep = "\t")

 #### Export Average Expression
 avg.oligos <- AverageExpression(dm, verbose = FALSE,group.by = c("group","clusters"))$RNA
 avg.oligos <- as_tibble(cbind(GeneLabel=rownames(avg.oligos),avg.oligos))
 avg.oligos <- filter(avg.oligos, GeneLabel %in% dm.markers$gene)
 Ranked <- data.frame(GeneLabel = arrange(dm.markers, avg_log2FC)$gene, stringsAsFactors = F)
 avg.oligos <- left_join(Ranked, avg.oligos)
 avg.oligos <- distinct(avg.oligos)
 write.table(avg.oligos,"RNA_AvgExpression_VarFeats_minExpessingCells10Pct_newSeurat.tsv",quote = F,row.names = F,col.names = T,sep = "\t")

 #### Get Abundance of expressing cells for every marker gene and group_cluster
 dm@meta.data$groupedClusters <- paste0(dm@meta.data$group,"_",dm@meta.data$clusters)
 Idents(dm) <- dm@meta.data$groupedClusters
 counts <- as.matrix(dm@assays$RNA@counts)
 counts <- data.frame(gene = row.names(counts),counts)
 counts <- filter(counts, gene %in% row.names(dm.markers))
 counts <- t(counts[-1])
 counts <- data.frame(cell = row.names(counts),groupedClusters = dm@meta.data$groupedClusters, counts)
 counts <- gather(counts,key = "gene",value = "count",-groupedClusters,-cell)
 counts$expressed <- counts$count > 0
 counts %>% group_by(groupedClusters,gene) %>% summarise(CellNumber = n(), AvgExpression = sum(count)/CellNumber, Pct1 = sum(expressed)/CellNumber) -> countStats
 countExport <- spread(select(countStats, -CellNumber, -Pct1), key = groupedClusters, value = AvgExpression)
 
 # Export mean expression
 write.table(countExport,"Counts_AvgExpression_VarFeats_newSeurat.tsv",quote = F,row.names = F,col.names = T,sep = "\t")

 # Export expressing cell abundance
 countPct <- spread(select(countStats, -CellNumber, -AvgExpression), key = groupedClusters, value = Pct1)
 write.table(countPct,"Counts_Pct1_VarFeats_newSeurat.tsv",quote = F,row.names = F,col.names = T,sep = "\t")




### Expression plots (Figure 4)
###############

### Plotting data import and preprocessing
 
### Import data for bar graphs
props <- read.csv("ClusterCell_Proportions_newSeurat.tsv",header = T,stringsAsFactors = F,sep="\t")
abundance <- read.csv("Counts_Pct1_AllFeats_newSeurat.tsv",header = T,stringsAsFactors = F,sep="\t")
abundance <- gather(abundance,key="group",value = "abundance",-gene)

# Add the number of cells in each cluster for weighting of the abundance of detected expression later
props$group_clusters <- paste0(props$group,"_",props$clusters)
groupAbundance <- left_join(abundance,select(props, group_clusters, cluster_cells),by = c('group' = 'group_clusters'))

abundance <- separate(abundance,col = "group",into = c("group","cluster"))
groupAbundance <- separate(groupAbundance,col = "group",into = c("group","cluster"))

# Infer the abundance of detected expression in the whole group independent of cluster
groupAbundance %>% group_by(gene,group) %>% summarise(abundance = sum((abundance * cluster_cells))/sum(cluster_cells),nclusters = n()) -> groupAbundance

# Relevel cluster names to be consistent with VlnPlots
abundance$cluster <- factor(abundance$cluster,levels = c("Neuronal","NPC","iOPC","iOL1","iOL2"))

### Plot function definitions

# Creates a barplot showing expressing cell abundance in a given cluster/sample
abundanceBarPlot <- function(feature){
  return(ggplot(filter(abundance,gene == feature))+
    geom_col(aes(x = cluster, y = 100*abundance ,fill = group, group = group, colour = "black"),position = position_dodge2(preserve = "single"))+
    theme_cowplot()+
    scale_fill_manual(values=c("gray","light blue","dark blue"),guide = F)+
    scale_colour_manual(values="black",guide=F)+
    ylab("Expressing Cells [%]")+
    ylim(0,100)+
    theme(axis.text.x = element_text(angle = 45,vjust = 1,hjust =1),plot.title.position = "panel")+
    xlab("Identity"))
}

# Creates a barplot showing overall expressing cell abundance in a given sample
groupAbundanceBarPlot <- function(feature){
  return(ggplot(filter(groupAbundance,gene == feature))+
           geom_col(aes(x = group, y = 100*abundance ,fill = group, group = group, colour = "black"),position = position_dodge2(preserve = "single",width = 0.2))+
           theme_cowplot()+
           scale_fill_manual(values=c("gray","light blue","dark blue"),guide = F)+
           scale_colour_manual(values="black",guide=F)+
           ylab("Expressing Cells [%]")+
           ylim(0,100)+
           theme(axis.text.x = element_text(angle = 45,vjust = 1,hjust =1),plot.title.position = "panel")+
           xlab("Identity"))
}

# Puts dimension, violin and bar plots together in a grid for export
featureGridPlot <- function(feature,label_size = 18){
  if(sum(abundance$gene == feature)){
    plot1 <- FeaturePlot(dm, repel = T, split.by = "group",sort.cell = T, features = feature,reduction = "umap",cols = c("light blue","dark red"),ncol = 1,pt.size = 1.5,label = T,label.size = 5,keep.scale = 'all')+theme(legend.position = "right")
    plot2 <- groupAbundanceBarPlot(feature)
    plot3 <- abundanceBarPlot(feature)
    plot4 <- VlnPlot(dm,feature,split.by = "group",cols = c("gray","light blue","dark blue"),pt.size = 0) + theme (legend.position = "top", plot.title = element_text(color = "white",size = 0), legend.spacing.x = unit(0.3,"line"))
    plot5 <- plot_grid(plot2,plot3,ncol=2,align = "hv",labels=c('B','C'),label_size = label_size,hjust = 1,rel_widths = c(1,2.5))
    plot6 <- plot_grid(plot5,plot4,nrow=2,labels=c('','D'),label_size = label_size,hjust = 1,rel_heights = c(1,1.2))
    plot7 <- plot_grid(plot1,plot6, align = "hv",axis = "l",labels = c('A'), label_size = label_size, rel_widths =c(3,1),rel_heights = c(2,1))
    title_theme <- ggdraw() + draw_label(feature, size = 20)
    plot8 <- plot_grid(title_theme, plot7, ncol = 1, rel_heights = c(0.2, 1))
    return(plot8)
  }
}

### Update sample name
groupAbundance$group[groupAbundance$group == "Pre"] <- "NI"
dm$group[dm$group == "Pre"] <- "NI Day+0"
dm$group[dm$group == "S"] <- "S Day+15"
dm$group[dm$group == "SON"] <- "SON Day+15"
props$group[props$group == "Pre"] <- "NI Day+0"
props$group[props$group == "S"] <- "S Day+15"
props$group[props$group == "SON"] <- "SON Day+15"

### Plotting loop
feature_list <- c("ASCL1", "BMI1", "HES1", "MSI1", "NES", "PAX6", "PROM1", "SOX1", "SOX2", "SOX3", "ST8SIA1", "PDGFRA", "CSPG4", "CD9", "CNP", "GPR17", "PTPRZ1", "MYT1", "OLIG1", "OLIG2", "ID2", "ID4", "SOX9", "SOX10", "CLDN11", "GALC", "MAG", "MAL", "MBP", "MOBP", "MOG", "MYRF", "NKX6-2", "NKX2-2", "PLP1", "SMARCA4", "SOX17", "TRF", "CD82", "ZFP191", "ZFP488", "ZFP536", "MKI67", "KCTD13", "ZNF536", "ZNF488", "MALAT1", "STMN2", "NEFM", "NEFL", "SNAP25", "NEUROD1", "NEUROD4", "SRY", "VIM", "SPARC", "SLC6A4", "GRIA1", "GRIA2", "GRIA3", "GRIA4", "GRIN1", "GRIN2A", "GRIN2B", "GRIN2C", "GRIN2D", "GRIN3A", "GRIN3B", "CNR1", "CNR2", "CACNA1C", "CACNA2D1", "CACNB2", "CDKN1C")
pdf("FeatureGridPlots_v05.pdf", width = 24, height = 8)
for(feature in feature_list){
  print(featureGridPlot(feature = feature))
}
dev.off()

# Reorder cluster levels
props$clusters <- factor(props$clusters,levels=c("Neuronal","NPC","iOPC","iOL1","iOL2"))

# Plot cluster proportions
ggplot(props)+
  geom_col(aes(x = clusters, y = percentage ,fill = group, group = group, colour = "black"),position = position_dodge2(preserve = "single"))+
  theme_cowplot()+
  scale_fill_manual(values=c("gray","light blue","dark blue"))+
  scale_colour_manual(values="black",guide=F)+
  ylab("Cell Type Abundance [%]")+
  theme(axis.text.x = element_text(angle = 45,vjust = 1,hjust =1))+
  ylim(0,100)+
  theme(legend.title = element_blank())+
  xlab("Identity")

###### Plot latent time
d <- dm@meta.data
d <- filter(d, clusters != "Neuronal")

### Get statistical metrics and test for constrast of interest (S vs. SON)
d %>% group_by(group) %>% summarise(mp = median(velocity_pseudotime), ml  = median(latent_time),mad=mad(latent_time),n=n())->md
wilcox.test(x = d$latent_time[d$group == "S Day+15"], y = d$latent_time[d$group == "SON Day+15"])

### Relevel group factor for plotting
d$group <- factor(d$group,levels = c("NI Day+0","S Day+15","SON Day+15"))

### Violin/Box plot
ggplot(d,aes(x=group,y=latent_time))+
  theme_cowplot()+
  ylab("Latent Time")+
  xlab("Sample")+
  ylim(0,1.25)+
  ggtitle("Latent time distribution")+
  geom_signif(comparisons = list(c("S Day+15","SON Day+15"),c("NI Day+0","SON Day+15"),c("S Day+15","NI Day+0")), y=c(1,1.1,1),map_signif_level = c("***" = 0.001, "**" = 0.01, "*" = 0.05,"#" = 0.1,"n.s." = 1.1), textsize = 8,size = 1.05,tip_length = 0.025)+
  geom_violin(aes(color = group,group = as.numeric(group)-0.2, fill = group,x = group, y = latent_time))+
  geom_boxplot(aes(group = group, x = as.numeric(group)+0.2, y = latent_time), color = "black",fill = "white",outlier.alpha = 0,size = 1.05,width = 0.1)+
  scale_color_manual(values=c("light grey","light blue","dark blue"))+
  scale_fill_manual(values=c("light grey","light blue","dark blue"))

### Split OPC cluster in two isolate them (Figure 6)

# Recluster
dm <- FindNeighbors(dm, dims = 1:pc)
dm <- FindClusters(dm, resolution = 0.28) ### Resolution increased from 0.23 to 0.28

# Test plot
print(DimPlot(dm,reduction = "umap", pt.size = 2,label = T))

# Update sample IDs
dm@meta.data$group <- factor(dm$group,levels = c("NI Day+0","S Day+15","SON Day+15"))

# Update cluster ids
new.cluster.ids <- c("iOPC1","iOPC2","NPC","Neuronal","iOL1","iOL2")
names(new.cluster.ids) <- levels(dm)
dm <- RenameIdents(dm, new.cluster.ids)
dm@meta.data$clusters2 <- unname(dm@active.ident)
dm@meta.data$group_clusters2 <- paste0(dm$group,"_",dm$clusters2)
newIDs <- dm$group_clusters2
names(newIDs) <- names(dm@active.ident)
Idents(dm) <- newIDs

do <- subset(dm, subset = clusters2 %in% c("iOPC1", "iOPC2"))
allFeats <- rownames(do@assays$SCT@data)
write.table(allFeats,"SCT_DGE_OPC_SplitClusters_BackgroundGeneset.tsv",quote = F,row.names = F,col.names = F,sep = "\t")


### Differential Expression Analysis OPC subpopulations
################

### Get full expressed gene list for OPC2 vs OPC1 comparisons
allFeats <- rownames(do@assays$SCT@data)
write.table(allFeats,"SCT_DGE_OLlineage_SplitClusters_BackgroundGeneset.tsv",quote = F,row.names = F,col.names = F,sep = "\t")



### S OPC2 vs. S OPC1

# Get DEGs
markers <- FindMarkers(dm_S,ident.1 = "S Day+15_iOPC2", ident.2 = "S Day+15_iOPC1", assay = "SCT",logfc.threshold = 0) # Setting the logfc threshold to 0 returns a full list of all expressed genes

# Annotate clusters compared
markers$Cluster1 <- "S Day+15_iOPC2"
markers$Cluster2 <- "S Day+15_iOPC1"

# Reformat dataset
markers$gene <- row.names(markers)
markers <- arrange(markers,desc(avg_log2FC),p_val_adj)

# Export
write.table(markers,"SCT_DGE_Sonly_OPC2vsOPC1_SplitClusters_unfiltered.tsv",quote = F,row.names = F,col.names = T,sep = "\t")

# Volcanoplot
VolcanoPlot(markers,log2FC_threshold = 1)+ggtitle("S Day+15 OPC2 vs. OPC1")

# Filter genelist for hits
markers <- filter(markers,p_val_adj < 0.05,avg_log2FC > 1 | avg_log2FC < -1 )

# Arrange list
markers <- arrange(markers,desc(avg_log2FC),p_val_adj)

# Export
write.table(markers,"SCT_DGE_Sonly_OPC2vsOPC1_SplitClusters_highLog2FC.tsv",quote = F,row.names = F,col.names = T,sep = "\t")



### SON OPC2 vs. S OPC2

markers <- FindMarkers(dm,ident.1 = "SON Day+15_OPC2", ident.2 = "S Day+15_OPC2", assay = "SCT",logfc.threshold = 0)
markers$Cluster1 <- "SON Day+15_OPC2"
markers$Cluster2 <- "S Day+15_OPC2"
markers$gene <- row.names(markers)
markers <- arrange(markers,desc(avg_log2FC),p_val_adj)
write.table(markers,"SCT_DGE_OPC2_SONvsS_SplitClusters_unfiltered_forGSEA.tsv",quote = F,row.names = F,col.names = T,sep = "\t")

VolcanoPlot(markers,log2FC_threshold = 1)+ggtitle("SON Day+15 OPC2 vs. S Day+15 OPC2")

markers <- filter(markers,p_val_adj < 0.05,abs(avg_log2FC) > 1)
markers <- arrange(markers,desc(avg_log2FC),p_val_adj)
write.table(markers,"SCT_DGE_OPC2_SONvsS_SplitClusters_highLog2FC.tsv",quote = F,row.names = F,col.names = T,sep = "\t")



### SON OPC2 vs S OPC1 (There are just a few SON OPC1 cells, producing a patchy and artefact prone expression profile)

markers <- FindMarkers(dm,ident.1 = "SON Day+15_OPC2", ident.2 = "S Day+15_OPC1", assay = "SCT",logfc.threshold = 0)
markers$Cluster1 <- "SON Day+15_OPC2"
markers$Cluster2 <- "S Day+15_OPC1"
markers$gene <- row.names(markers)
markers <- arrange(markers,desc(avg_log2FC),p_val_adj)
write.table(markers,"SCT_DGE_SON_OPC2vsS_OPC1_SplitClusters_unfiltered_forGSEA.tsv",quote = F,row.names = F,col.names = T,sep = "\t")

VolcanoPlot(markers,log2FC_threshold = 1,adjust_p_val = T)+ggtitle("SON Day+15 OPC2 vs. S Day+15 OPC1")

markers <- filter(markers,p_val_adj < 0.05,abs(avg_log2FC) > 1)
markers <- arrange(markers,desc(avg_log2FC),p_val_adj)
write.table(markers,"SCT_DGE_SON_OPC2vsS_OPC1_SplitClusters_highLog2FC.tsv",quote = F,row.names = F,col.names = T,sep = "\t")