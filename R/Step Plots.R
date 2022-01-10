#############   Construct Expression Ranking Keppler-Müller-Plots #############
# Version:      1
# Authour:      Marius Stephan
# Created:      22.09.2021
# Last edited:  22.09.2021

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

# Define local functions
getRanking <- function(d,desc=F){
  d <- data.frame(d=d,p=c(1:length(d)))
  d <- arrange(d,d)
  d$r <- c(1:length(d$d))
  d <- arrange(d,p)
  return(d$r)
}

graph.scatter <- function (data.df,                                       # Dataset as data.frame or tibble
                           grouping.var = "Group",                        # Name of independent variable
                           x,                                             # Dependent variable mapped onto x axis
                           y,                                             # Dependent variable mapped onto y axis
                           
                           omit.NA = TRUE,                                # There is virtually no reason to turn this off. NA values are of no use for plotting whatsoever
                           
                           legend.pos=c("right"),                         # Coordinates of the legend c(x,y). Can be customized (see ggplot2 Doc for details)
                           text.size=20,                                  # Text size of any labels. Default is set to make the graph readable at low resolution
                           x.label = x,                                   # Label of x-axis and legend (should represent all groups)
                           y.label = y,                                   # Label of the y-axis. Default is the dep. variable name.
                           graph.title = NULL,                            # Title of the graph for later identification. Default set to sth. annoying to force people to think about it.
                           custom.colours = NULL,                         # A vector of colours (either named, hex, cmyk, rgb or hcl, you name it) assigned to each group in plotting order. Default is c("dark grey","dark blue","grey","light blue")
                           alpha = 1,                                     # Opacity of datapoints
                           
                           # Spearman Correlation Analysis
                           correlation = FALSE,                           # Plot regression line of Spearman correlation analysis? Default false.
                           corr.label = TRUE,                             # Plot R² and p-value (Strongly recommended; the regression line alone yields limited information)
                           corr.se = FALSE){                              # Display SE-Ribbon in back ground (recommended for ONE regression ONLY)
  
  ### XY Scatter Plot (w/ Spearman correlation)
  ### ________________________________________________________________________________________________
  
  #---------------------------------------------------
  ## Rename variables
  #___________________________________________________
  #---------------------------------------------------
  
  data.df <- rename(data.df, group = first(grouping.var), x = first(x), y = first(y))
  
  
  #---------------------------------------------------
  ## Drop NAs
  #___________________________________________________
  #---------------------------------------------------
  
  if (omit.NA == T) {
    data.df <- data.df[!is.na(data.df$x && !is.na(data.df$y)),]
  }
  
  
  #---------------------------------------------------                   
  ## Group Data
  #___________________________________________________
  #---------------------------------------------------
  
  
  #---------------------------------------------------
  ## Get group names
  #---------------------------------------------------
  
  group.names <- sort(unique(pull(data.df,group)))
  
  #---------------------------------------------------
  ## Join data.df and grouping
  #---------------------------------------------------
  
  group.tbl <- c()
  group.tbl$group <- sort(unique(pull(data.df,group)))
  group.tbl$order <- c(1:length(group.tbl$group))
  
  group.tbl <- as_tibble(group.tbl)
  data.df <- inner_join(data.df,group.tbl, by = "group")
  data.df$order <- as.character(data.df$order)
  data.df <- data.df[order(data.df$order),]
  
  
  #---------------------------------------------------
  ## Calculate Spearman Correlation
  #___________________________________________________
  #---------------------------------------------------
  
  data.df %>% 
    group_by(group, order) %>%
    summarise(r = c(unname(cor.test(x,y)[[4]]),0),
              r2 = c(round(r[1]**2,3),0),
              p = c(round(cor.test(x,y)[[3]],3),0),
              label = c( sprintf("~R^{2} == '%.2f'",r2[1]),sprintf("~p == '%.3f'",p[1]))) -> cor.df
  
  #--------------------------------------------------- 
  ## Determine other plotting modalities
  #___________________________________________________
  #---------------------------------------------------  
  
  #---------------------------------------------------
  ## Check parameters for custom title
  #---------------------------------------------------
  if (is.null(graph.title)) {
    graph.title = ""  # Blanks ggtitle, if no custom title has been given.
  }
  
  #---------------------------------------------------
  ## Determine Axis Limits
  #---------------------------------------------------
  
  # X-Limits
  
  xlimits <- c( min(data.df$x) - 0.1*max(data.df$x), 1.1*max(data.df$x))
  
  # Y-limits
  
  ylimits <- c( min(data.df$y) - 0.1*max(data.df$y), 1.2*max(data.df$y))
  
  
  #---------------------------------------------------
  ## Determine plotting colours
  #---------------------------------------------------
  group.number <- length(group.names)
  plot.colours <- graph.helper.colourFinder(group.number, custom.colours)
  
  plot <- ggplot(data.df,aes(x=x, y=y,colour=order))+
    
    #-------------------------------------------------
  ## Format axes, labels, background and legend
  #-------------------------------------------------
  
  # Add title
  
  ggtitle(graph.title)+
    
    # Customize axis labels
    
    labs(colour="black",x=x.label,y=y.label)+
    
    # Make all unwanted elements blank
    
    theme_bw()+
    theme(panel.grid.major = element_blank(),
          panel.border = element_blank(),
          panel.grid.minor = element_blank(),
          legend.key = element_blank(),
          strip.background = element_blank(),
          strip.text = element_blank(),
          
          # Customize axes
          
          axis.line = element_line(colour = "black", size = 1.1),
          axis.ticks = element_line(size = 1.1),
          axis.text = element_text(colour = "black"),
          text= element_text(size = text.size, colour = "black"),
          
          # Customize legend
          
          legend.position = legend.pos,
          legend.key.height = unit(3,"line"),
          legend.text = element_text(size = text.size))+
    
    scale_colour_manual(labels = group.names, name = "", values = plot.colours)+
    scale_fill_manual(values = plot.colours,guide = "none")+
    scale_x_continuous(limits=xlimits)+
    scale_y_continuous(limits=ylimits)+
    geom_point(alpha=alpha,size=1.5)
  
  if (correlation) {
    plot <- plot + geom_smooth(aes(fill=order),method=lm, se=corr.se, fullrange=T,formula = y ~ x,size=1.1,show.legend = F)
    if (corr.label) {
      plot <- plot + geom_text(aes(x=xlimits[2],y=ylimits[2],group=order, label=cor.df$label, vjust=sort(c(1.2*seq(1,group.number,1)-1.1, 1.2*seq(1,group.number,1)-0.8)), hjust = rep(c(2,1),group.number)), parse = T, data = cor.df,size=text.size/3,show.legend = FALSE ) +
        coord_cartesian(clip = 'off')
    }
  }
  
  return(plot)
}


# Load user-defined functions
source("C:/Users/Marius/OneDrive - LMU Uni München/01.Stephan/RStudio Projects/My Functions/graph.functions v21.R")

# Load data
dm <- readRDS("OPCmerge_v09_OPCsubclusters.rds")

# Get the meta data table
dm@meta.data$groupedClusters <- paste0(dm@meta.data$group,"_",dm@meta.data$clusters)

# Define marker geneset for marker load
markerSet <- c("CNP", "CD9", "CLDN11", "GALC", "PLP1", "MAG", "MAL", "MBP", "MOBP", "MYRF")

# Calculate marker load (more reliable than single features)
t <- as.data.frame(t(dm@assays$SCT@data[row.names(dm@assays$SCT@data) %in% markerSet,]))
t$normMarkerLoad <- rowSums(t)
t$cell <- row.names(t)
t <- arrange(t,cell)

# Add to meta.data table
dm@meta.data$cell <- row.names(dm@meta.data)
dm@meta.data <- arrange(dm@meta.data,cell)
dm@meta.data <- cbind(dm@meta.data,select(t,normMarkerLoad))


##########
### Add ranking to each sample
##########

addRanks2Samples <- function(d){
  # Extract relevant data from meta.data table
  S <- select(filter(d@meta.data, group == "S Day+15"),cell, construct_normExpression, normMarkerLoad)
  SON <- select(filter(d@meta.data, group == "SON Day+15"),cell, construct_normExpression, normMarkerLoad)
  
  # Rank cells by construct expression and marker load
  S$normRank <- getRanking(S$construct_normExpression)
  S$normMarkerRank <- getRanking(S$normMarkerLoad)
  S$normCumulative <- S$normRank/max(S$normRank)
  S$normMarkerCumulative <- S$normMarkerRank/max(S$normMarkerRank)
  
  SON$normRank <- getRanking(SON$construct_normExpression)
  SON$normMarkerRank <- getRanking(SON$normMarkerLoad)
  SON$normCumulative <- SON$normRank/max(SON$normRank)
  SON$normMarkerCumulative <- SON$normMarkerRank/max(SON$normMarkerRank)
  
  # Add ranks to meta.data
  t <- select(arrange(select(rbind(S,SON),cell,normRank,normMarkerRank,normCumulative,normMarkerCumulative),cell),-cell)
  d@meta.data <- arrange(d@meta.data,cell)
  d@meta.data <- cbind(d@meta.data,t)
  
  return(d)
}

###########
### iOL+ iOPCs
iOall <- subset(dm, clusters %in% c("OPC","OL1","OL2") & group != "NI Day+0")
iOall <- addRanks2Samples(iOall)


pdf("Step Functions v06.pdf",height=5,width=7)

# Kolmogorov-Smirnow test for distribution differences
test <- ks.test(x = subset(iOall,group == "S Day+15")@meta.data$construct_normExpression,y = subset(iOall,group == "SON Day+15")@meta.data$construct_normExpression)

# Wilcoxon test alternative
# wilcox.test(x = subset(iOall,group == "S Day+15")@meta.data$construct_normExpression,y = subset(iOall,group == "SON Day+15")@meta.data$construct_normExpression)

y.max <- max(iOall@meta.data$construct_normExpression)
x.max <- max(iOall@meta.data$normCumulative)
ggplot(arrange(iOall@meta.data,normCumulative))+ 
  ggtitle("Construct Expression iOalls")+
  geom_line(aes(y=construct_normExpression,x=normCumulative*100,group=group,colour=group),size=2)+
  geom_text(data=data.frame(x=c(1)),mapping=aes(x=50,y=1.2*y.max),label=sprintf("D = %.2f   p = %.3f",test[[1]],test[[2]]),size =6)+
  scale_color_manual(values = c("light blue", "dark blue"))+
  theme(element_text(size=28))+
  #ylim(0,120)+
  ylab("Construct Expression")+
  xlab("Cumulaticve Cell Count [%]")+
  theme_cowplot()+
  theme(axis.title = element_text(size=18), axis.text = element_text(size = 18), axis.ticks = element_line(size = 1.4),axis.line = element_line(size = 1.4),
        legend.text = element_text(size=18),legend.title = element_blank(),legend.position = "right")

test <- ks.test(x = subset(iOall,group == "S Day+15")@meta.data$normMarkerLoad,y = subset(iOall,group == "SON Day+15")@meta.data$normMarkerLoad)
y.max <- max(iOall@meta.data$normMarkerLoad)
x.max <- max(iOall@meta.data$normMarkerCumulative)

ggplot(arrange(iOall@meta.data,normMarkerCumulative))+ 
  ggtitle("Marker Load iOL Lineage")+
  geom_step(aes(y=normMarkerLoad,x=normMarkerCumulative*100,group=group,colour=group),size=2,)+
  geom_text(data=data.frame(x=c(1)),mapping=aes(x=50,y=1.2*y.max),label=sprintf("D = %.2f   p = %.3f",test[[1]],test[[2]]),size =6)+
  scale_color_manual(values = c("light blue", "dark blue"))+
  ylab("Marker Load")+
  xlab("Cumulaticve Cell Count [%]")+
  theme_cowplot()+
  theme(axis.title = element_text(size=18), axis.text = element_text(size = 18), axis.ticks = element_line(size = 1.4),axis.line = element_line(size = 1.4),
        legend.text = element_text(size=18),legend.title = element_blank(),legend.position = "right")

graph.scatter(iOall@meta.data,x.label = "Marker Load", y.label = "Construct Expression",graph.title = "iOall Constr. Expr. vs Marker Load (norm.)",text.size = 18, grouping.var = "group", x = "normMarkerLoad", y = "construct_normExpression",correlation = T,custom.colours = c("light blue", "dark blue"),corr.se = T)


dev.off()

