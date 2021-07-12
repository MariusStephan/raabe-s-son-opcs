#############   DGE Venn Diagrams  #############
# Version:      01
# Authour:      Marius Stephan
# Created:      18.05.2021
# Last edited:  18.05.2021

# Clear Workspace
rm(list=ls())

# Load Packages needed
library(tidyverse)
library(ggvenn)
library(ggVennDiagram)

# Import data
up <- list( S_OPC2vsS_OPC1=filter(read.delim("SCT_DGE_Sonly_OPC2vsOPC1_SplitClusters_unfiltered_forGSEA.tsv"),p_val_adj < 0.05, avg_log2FC > 1)$gene,
            SON_OPC2vsS_OPC1=filter(read.delim("SCT_DGE_SON_OPC2vsS_OPC1_SplitClusters_unfiltered_forGSEA.tsv"),p_val_adj < 0.05, avg_log2FC > 1)$gene
          )
down <- list( S_OPC2vsS_OPC1=filter(read.delim("SCT_DGE_Sonly_OPC2vsOPC1_SplitClusters_unfiltered_forGSEA.tsv"),p_val_adj < 0.05, avg_log2FC < -1)$gene,
              SON_OPC2vsS_OPC1=filter(read.delim("SCT_DGE_SON_OPC2vsS_OPC1_SplitClusters_unfiltered_forGSEA.tsv"),p_val_adj < 0.05, avg_log2FC < -1)$gene
            )

pdf("VennDiagrams_OPC2vsOPC1.pdf",height = 4, width = 5)
### Plot Venn diagrams
ggVennDiagram(up,label_alpha = 0,label_size = 6)+ggtitle("Upregulate")+scale_fill_gradient(low = "lightpink",high = " dark red")+scale_color_manual(values = c("white","white"))
ggVennDiagram(down,label_alpha = 0,label_size = 6)+ggtitle("Downregulated")+scale_fill_gradient(low = "light blue",high = " dark blue")+scale_color_manual(values = c("white","white"))
dev.off()

up <- list( OL2vsOPC1=filter(read.delim("SCT_DGE_OL2vsOPC1_SplitClusters_unfiltered_forGSEA.tsv"),p_val_adj < 0.05, avg_log2FC > 1)$gene,
            OL2vsOPC2=filter(read.delim("SCT_DGE_OL2vsOPC2_SplitClusters_unfiltered_forGSEA.tsv"),p_val_adj < 0.05, avg_log2FC > 1)$gene
)
down <- list( OL2vsOPC1=filter(read.delim("SCT_DGE_OL2vsOPC1_SplitClusters_unfiltered_forGSEA.tsv"),p_val_adj < 0.05, avg_log2FC < -1)$gene,
              OL2vsOPC2=filter(read.delim("SCT_DGE_OL2vsOPC2_SplitClusters_unfiltered_forGSEA.tsv"),p_val_adj < 0.05, avg_log2FC < -1)$gene
)

pdf("VennDiagrams_OL2vsOPCs.pdf",height = 4, width = 5)
ggVennDiagram(up,label_alpha = 0,label_size = 6)+ggtitle("Upregulated")+scale_fill_gradient(low = "lightpink",high = " dark red")+scale_color_manual(values = c("white","white"))
ggVennDiagram(down,label_alpha = 0,label_size = 6)+ggtitle("Downregulated")+scale_fill_gradient(low = "light blue",high = " dark blue")+scale_color_manual(values = c("white","white"))
dev.off()


### Lowered threshold
up <- list( S_OPC2vsS_OPC1=filter(read.delim("SCT_DGE_Sonly_OPC2vsOPC1_SplitClusters_unfiltered_forGSEA.tsv"),p_val_adj < 0.05, avg_log2FC > 0.25)$gene,
            SON_OPC2vsS_OPC1=filter(read.delim("SCT_DGE_SON_OPC2vsS_OPC1_SplitClusters_unfiltered_forGSEA.tsv"),p_val_adj < 0.05, avg_log2FC > 0.25)$gene
)
down <- list( S_OPC2vsS_OPC1=filter(read.delim("SCT_DGE_Sonly_OPC2vsOPC1_SplitClusters_unfiltered_forGSEA.tsv"),p_val_adj < 0.05, avg_log2FC < -0.25)$gene,
              SON_OPC2vsS_OPC1=filter(read.delim("SCT_DGE_SON_OPC2vsS_OPC1_SplitClusters_unfiltered_forGSEA.tsv"),p_val_adj < 0.05, avg_log2FC < -0.25)$gene
)

pdf("VennDiagrams_OPC2vsOPC1_lowFCthreshold.pdf",height = 4, width = 5)
ggVennDiagram(up,label_alpha = 0,label_size = 6)+ggtitle("Upregulate")+scale_fill_gradient(low = "lightpink",high = " dark red")+scale_color_manual(values = c("white","white"))
ggVennDiagram(down,label_alpha = 0,label_size = 6)+ggtitle("Downregulated")+scale_fill_gradient(low = "light blue",high = " dark blue")+scale_color_manual(values = c("white","white"))
dev.off()

up <- list( OL2vsOPC1=filter(read.delim("SCT_DGE_OL2vsOPC1_SplitClusters_unfiltered_forGSEA.tsv"),p_val_adj < 0.05, avg_log2FC > 0.25)$gene,
            OL2vsOPC2=filter(read.delim("SCT_DGE_OL2vsOPC2_SplitClusters_unfiltered_forGSEA.tsv"),p_val_adj < 0.05, avg_log2FC > 0.25)$gene
)
down <- list( OL2vsOPC1=filter(read.delim("SCT_DGE_OL2vsOPC1_SplitClusters_unfiltered_forGSEA.tsv"),p_val_adj < 0.05, avg_log2FC < -0.25)$gene,
              OL2vsOPC2=filter(read.delim("SCT_DGE_OL2vsOPC2_SplitClusters_unfiltered_forGSEA.tsv"),p_val_adj < 0.05, avg_log2FC < -0.25)$gene
)

pdf("VennDiagrams_OL2vsOPCs_lowFCthreshold.pdf",height = 4, width = 5)
ggVennDiagram(up,label_alpha = 0,label_size = 6)+ggtitle("Upregulated")+scale_fill_gradient(low = "lightpink",high = " dark red")+scale_color_manual(values = c("white","white"))
ggVennDiagram(down,label_alpha = 0,label_size = 6)+ggtitle("Downregulated")+scale_fill_gradient(low = "light blue",high = " dark blue")+scale_color_manual(values = c("white","white"))
dev.off()


