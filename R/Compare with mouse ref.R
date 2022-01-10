#############   Correlate scRNA-Seq data  with bulk RNA-Seq data based on Nat Biotech Paper #############
# Version:      02
# Authour:      Marius Stephan
# Created:      26.01.2021
# Last edited:  24.11.2021

# Clear Workspace
rm(list=ls())

# Load Packages needed
library(tidyverse)
library(xlsx)
library(DESeq2)
library(limma)
library(RColorBrewer)
library(pheatmap)
library(cowplot)
library(ggrepel)


### Site of the reference dataset used
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE52564

### Download link for human/mouse homologues genes list
# https://www.informatics.jax.org/downloads/reports/HOM_MouseHumanSequence.rpt

### Local Defined Functions
importXLS <- function(file){
  # Extract sample name
  n <- str_remove(str_split(file,"_",simplify = T)[2],".xls")
  d <- read.xlsx(file = file, sheetIndex = 1, header = T)
  colnames(d) <- c("Symbol", n)
  return(d)
}

importAllXLS <- function(path){
  cdir <- getwd()
  setwd(path)
  temp = list.files(pattern="*.xls")
  myfiles = lapply(temp, importXLS)
  setwd(cdir)
  
  return(myfiles)
}


### Create human/mouse ref table
r <- read.table("HOM_MouseHumanSequence.rpt", sep = "\t", header = T,stringsAsFactors = F)
r$Common.Organism.Name[r$Common.Organism.Name == "mouse, laboratory"] <- "mouse"
r %>% group_by(HomoloGene.ID) %>% summarise(Human = sum(Common.Organism.Name == "human"), Mouse = sum(Common.Organism.Name == "mouse")) -> r2
r2 <- filter(r2,Human == 1 & Mouse == 1)
r3 <- filter(r, HomoloGene.ID %in% r2$HomoloGene.ID)
r3 %>% group_by(HomoloGene.ID) %>% summarise(human = Symbol[Common.Organism.Name == "human"],mouse = Symbol[Common.Organism.Name == "mouse"]) -> h
rm("r","r2","r3")

### Import ref data
l <- importAllXLS("PATH_TO_DOWNLOADED_REF_DATA/RefData/")

### Join samples in one data.frame
d <- l[[1]]
print(paste0("1 of ", length(l), " files processed!"))
for (i in 2: length(l)){
  d <- left_join(d,l[[i]], by = "Symbol")
  print(paste0(i, " of ", length(l), " files processed!"))
}

### Save reference data to text file
write.table(d,"RefData.tsv", quote = F, sep="\t", col.names = T, row.names = F)

### Read reference data from file
r <- read.csv("RefData.tsv",header = T, sep = "\t", stringsAsFactors = F)
d <- read.csv("Counts_AvgExpression_VarFeats_newSeurat.tsv",header = T, sep = "\t", stringsAsFactors = F)

### Filter celltypes
r <- select(r, !starts_with(c("Microglia","Astrocyte","WC","Endothelial")))

d <- select(d, -NI.Day.0_OPC) # Based on 3 cells only, giving a very patchy expression profile
#d <- select(d, !ends_with("Neuronal"))
colnames(d) <- c("GeneLabel",colnames(d)[-1])
colnames(d) <- str_replace_all(colnames(d),"_O","_iO")

### Filter genes for homologues and intersect
i <- intersect(d$GeneLabel,h$human)
d <- filter(d, GeneLabel %in% i)

i <- h$mouse[h$human %in% i]
r <- filter(r, Symbol %in% i)

r <- left_join(r,h[2:3], by = c("Symbol" = "mouse"))
r$GeneLabel <- r$human
r <- select(r, -human, -Symbol)
  
### Create a DESeq2 object
#############

### Reformat datasets
i <- intersect(r$GeneLabel,d$GeneLabel)
d <- filter(d,GeneLabel %in% i)
r <- filter(r,GeneLabel %in% i)

d2 <- inner_join(d,r, by = "GeneLabel")
row.names(d2) <- d2$GeneLabel
d2 <- as.data.frame(select(d2,-GeneLabel))

### Provide Column Data
colData <- data.frame(condition = c(colnames(d)[-1], str_remove(colnames(r)[-length(r)],"[1-9]")), experiment = c(rep("1",length(d)-1),rep("2",length(r)-1)), type = rep("single-read",length(d2)))

### Adjust count matrix
d3 <- round(d2*1e4)
# Note: AvgExpression() computes the geometric mean of UMI counts, which doesn't produce integers.
# DESeq2 expects raw counts, so we needed to round the averaged counts and increased precision by for 4 orders of magnitude.
# Latter doesn't change the result notably.

d3 <- as.data.frame(d3)
for (i in 1:length(d3)){ ### This loop takes care of a bug, we've experienced.
  d3[i] <- as.integer(unlist(d3[i]))
}

### DESeq2 object from count matrix
dds <- DESeqDataSetFromMatrix(countData = d3, colData = colData, design = ~ condition)
vst <- varianceStabilizingTransformation(dds)
dat <- assay(vst)
dat2 <- removeBatchEffect(dat,colData$experiment) # Sample groups are linearly additive, so this won't work well if at all. Doesn't hurt either. The alternative option is to unblind DESeq2's modelling algorhithm in DESeqDataSetFromMatrix(). Thankfully, that step wasn't necessary in this pipeline.
assay(vst) <- dat2

### Check for major batch effects between experiments
plotPCA(vst)

### Compute gene expression variance across samples
rv <- rowVars(assay(vst))

### Select the 50 most variable genes for plotting
select <- order(rv, decreasing = TRUE)[seq_len(min(50, length(rv)))]

### Get the DESeqTransform object with the top 50 most variable genes
vst50 <- vst[select,]
plotPCA(vst50)

### Select the 2000 most variable genes to go into the comparison.
# This is done to improve SNR, but not entirely necessary as all features were preselected to avoid excessive noise. Check the methods section for details.
select <- order(rv, decreasing = TRUE)[seq_len(min(2000, length(rv)))]
vst2000 <- vst[select,]

### Plot PCA dimension plot with selected features only (see Figure 5).
update_geom_defaults("point",list(size=5))

pcaData <- plotPCA(vst2000, intgroup=c("condition", "type"), returnData = T)
percentVar <- round(100 * attr(pcaData, "percentVar"))
pcaData$name <- str_remove(str_remove(str_remove(pcaData$name,".Day."),"15"),"0")
ggplot(pcaData, aes(PC1, PC2, colour = condition)) +
  geom_point(size=3, show.legend = F) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()+
  geom_text_repel(aes(label = name),show.legend = F)+
  theme_cowplot()+
  scale_color_manual(values= c(rgb(0,0,0),rgb(0.7,0.7,0.7),rgb(0.3,0.3,0.3),hsv(0.16,0.6,0.6),hsv(0.25,0.4,0.8),rgb(0.5,0.5,0.5),hsv(0.6,0.6,0.7),hsv(0.6,1,0.6),hsv(0.6,0.4,0.8),hsv(0.16,0.6,0.6),hsv(0.6,0.6,0.7),hsv(0.6,1,0.6),hsv(0.6,0.4,0.8),hsv(0.16,0.6,0.6)))

### Compute manhattan distances
distances_samples <- dist(t(assay(vst2000)), method = "manhattan")
distances_genes <- dist(assay(vst50), method = "manhattan")

### Reformat input data
dm <- as.matrix(log2(assay(vst50)))
dms <- as.matrix(distances_samples)
rownames(dm) <- rownames(assay(vst50))
colnames(dm) <- colnames(assay(vst50))
rownames(dms) <- vst$condition
colnames(dms) <- NULL

### Plot sample distances
colors <- colorRampPalette(colors = rev(c(rgb(1,1,1),rgb(0.9,0.9,1),rgb(0.2,0.2,0.6),rgb(0,0,0.4))),bias=0.7,) (400)
colors <- c(colors, rep("#FFFFFF",100))
pheatmap(dms,clustering_distance_rows = distances_samples,clustering_distance_cols =distances_samples,col=colors,lwd=1.2)

### Plot expression of the top 50 most variable genes
colors <- colorRampPalette(colors = c(rgb(1,1,1),rgb(0.9,0.9,1),rgb(0.2,0.2,0.6),rgb(0,0,0.4)),bias=1,) (400)
colors <- c( rep("#FFFFFF",100),colors)
pheatmap(dm,clustering_distance_rows = distances_genes,clustering_distance_cols = distances_samples,col=colors,lwd=1.2)

### Spearman correlation
c <- as.dist(cor(assay(vst2000),method = "spearman"))
cm <- as.matrix(c)
colors <- colorRampPalette(colors = c(rgb(1,1,1),rgb(0.9,0.9,1),rgb(0,0,0.4)),bias=1,) (500)
colors <- c( rep("#FFFFFF",200),colors)
pheatmap(cm,clustering_distance_rows = distances_samples,clustering_distance_cols = distances_samples,col=colors,lwd=1.2)
