## Heatmaps of abundant pathways at each stage

v <- read.csv("abund_path_veg.csv",header=T,row.names=1)
dim(v)
path <- read.csv("pathway.csv",header=T,row.names=1) 
data <- log2(path + 1)
dim(data)
keep <- rownames(v)
meta <- read.csv("metadata.csv",header=T)
dim(meta)
meta$Sample = colnames(data)
rownames(meta) = colnames(data)
library(tidyverse)
library(colorRamp2)
library(ComplexHeatmap)
library(RColorBrewer)

x <- data[rownames(data), meta %>% 
            filter(Developmental_stage=='Vegetative') %>% pull(Sample)]
dim(x)
data_subset <- x[rownames(x) %in% keep, ]
dim(data_subset)

heat <- t(scale(t(data_subset)))

myBreaks <- seq(-3, 3, length.out = 100)
myCol <- colorRampPalette(brewer.pal(9, "YlGnBu"))(100)

h1 <- Heatmap(heat,                                                 
              col=colorRamp2(myBreaks, myCol),
              show_heatmap_legend = FALSE,
              cluster_rows=TRUE,
              show_row_dend=FALSE,
              row_title="Pathway",
              row_title_side="left",
              row_title_gp=gpar(fontsize=14, fontface="bold"),
              row_title_rot=90,                               
              show_row_names=TRUE,
              row_names_side="left",
              row_names_gp=gpar(fontsize=10),                            
              cluster_columns=TRUE,                           
              show_column_dend=FALSE,                         
              column_title="Vegetative",                      
              column_title_side="top",
              column_title_gp=gpar(fontsize=16, fontface="bold"),
              column_title_rot=0,                                
              show_column_names=FALSE,                           
              width = unit(9, "cm"))



r <- read.csv("abund_path_rep.csv",header=T,row.names=1)
keep <- rownames(r)
x <- data[rownames(data), meta %>% 
            filter(Developmental_stage=='Reproductive') %>% pull(Sample)]

dim(x)
data_subset <- x[rownames(x) %in% keep, ]
dim(data_subset) 
heat <- t(scale(t(data_subset)))

h2 <- Heatmap(heat,                                                 
              col=colorRamp2(myBreaks, myCol),
              show_heatmap_legend = FALSE,
              cluster_rows=TRUE,
              show_row_dend=FALSE,
              row_title= NA,                               
              show_row_names=TRUE,
              row_names_side="left",
              row_names_gp=gpar(fontsize=10),                            
              cluster_columns=TRUE,                           
              show_column_dend=FALSE,                         
              column_title="Reproductive",                      
              column_title_side="top",
              column_title_gp=gpar(fontsize=16, fontface="bold"),
              column_title_rot=0,                                
              show_column_names=FALSE,                           
              width = unit(9, "cm"))



m <- read.csv("abund_path_mat.csv",header=T,row.names=1)
keep <- rownames(m)
x <- data[rownames(data), meta %>% 
            filter(Developmental_stage=='Maturation') %>% pull(Sample)]

dim(x)
data_subset <- x[rownames(x) %in% keep, ]
dim(data_subset) 
heat <- t(scale(t(data_subset)))

h3 <- Heatmap(heat,                                                 
              col=colorRamp2(myBreaks, myCol),
              name="Abundance Z Score",
              heatmap_legend_param=list(
                color_bar="continuous",
                legend_direction="vertical",
                legend_width=unit(5,"cm"),
                title_position="leftcenter-rot",
                title_gp=gpar(fontsize=15, fontface="bold")),
              cluster_rows=TRUE,
              show_row_dend=FALSE,
              row_title= NA,                               
              show_row_names=TRUE,
              row_names_side="left",
              row_names_gp=gpar(fontsize=10),                            
              cluster_columns=TRUE,                           
              show_column_dend=FALSE,                         
              column_title="Maturation",                      
              column_title_side="top",
              column_title_gp=gpar(fontsize=16, fontface="bold"),
              column_title_rot=0,                                
              show_column_names=FALSE,                           
              width = unit(9, "cm"))


library(ggplotify)
library(grid)
library(cowplot)
library(colorspace) 
h1 <- as.ggplot(h1)
h2 <- as.ggplot(h2)
h3 <- as.ggplot(h3)     
plot_grid(h1,h2,h3,ncol=3, labels=LETTERS[1:3], rel_widths = c(1, 1))    




tiff("abundant_pathways.tiff",height=10400,width=12600,res=600)
plot_grid(h1,h2,h3,ncol=3, labels=LETTERS[1:3], rel_widths = c(1, 1)) 
dev.off()   
