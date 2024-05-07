## Heatmaps showing enriched bacterial genera at each developmental stage

library(ComplexHeatmap)
library(RColorBrewer)
library(colorRamp2)
count <- read.csv("relative_abundance_bacteria.csv",header=T,row.names=1)
meta <- read.csv("metadata.csv",header=T,row.names=1)

rel <- read.csv("enriched_bacteria_veg_stage.csv",header=T,row.names=1)
keep <- rownames(rel)
data_subset <- count[rownames(count) %in% keep, ]
developmental_stage <- "Vegetative"  # Specify the developmental stage of interest
meta_stage <- subset(meta, Developmental_stage == developmental_stage)
count_stage <- data_subset[, colnames(data_subset) %in% rownames(meta_stage)]               
dim(count_stage)

myCol <- colorRampPalette(brewer.pal(7, "RdPu"))(35)
myBreaks <- seq(-3, 3, length.out = 35)
f <- t(scale(t(count_stage)))
h1 <- Heatmap(f,                                                 
              col=colorRamp2(myBreaks, myCol),
              show_heatmap_legend = FALSE,
              cluster_rows=TRUE,
              show_row_dend=FALSE,
              row_title="Bacteria",
              row_title_side="left",
              row_title_gp=gpar(fontsize=22, fontface="bold"),
              row_title_rot=90,                               
              show_row_names=TRUE, row_names_gp = gpar(fontsize = 20, fontface="bold"),                           
              cluster_columns=TRUE,                           
              show_column_dend=FALSE,                         
              column_title="Vegetative",                      
              column_title_side="top",
              column_title_gp=gpar(fontsize=22, fontface="bold"),
              column_title_rot=0,                                
              show_column_names=FALSE,                           
              width = unit(9, "cm"))

rel <- read.csv("enriched_bacteria_rep_stage.csv",header=T,row.names=1)
keep <- rownames(rel)
data_subset <- count[rownames(count) %in% keep, ]
developmental_stage <- "Reproductive"  # Specify the developmental stage of interest
meta_stage <- subset(meta, Developmental_stage == developmental_stage)
count_stage <- data_subset[, colnames(data_subset) %in% rownames(meta_stage)]               
dim(count_stage)
f <- t(scale(t(count_stage)))
h2 <- Heatmap(f,                                                 
              col=colorRamp2(myBreaks, myCol),
              show_heatmap_legend = FALSE,
              cluster_rows=TRUE,
              show_row_dend=FALSE,
              row_title=NA,        
              
              show_row_names=TRUE, row_names_gp = gpar(fontsize = 20, fontface="bold"),                           
              cluster_columns=TRUE,                           
              show_column_dend=FALSE,                         
              column_title="Reproductive",                      
              column_title_side="top",                          
              column_title_gp=gpar(fontsize=22, fontface="bold"),
              column_title_rot=0,                                
              show_column_names=FALSE,                           
              width = unit(9, "cm"))      




rel <- read.csv("enriched_bacteria_mat_stage.csv",header=T,row.names=1)
keep <- rownames(rel)
data_subset <- count[rownames(count) %in% keep, ]
developmental_stage <- "Maturation"  # Specify the developmental stage of interest
meta_stage <- subset(meta, Developmental_stage == developmental_stage)
count_stage <- data_subset[, colnames(data_subset) %in% rownames(meta_stage)]               
dim(count_stage)
f <- t(scale(t(count_stage)))
h3 <- Heatmap(f,                                                 
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
              row_title=NA,        
              
              show_row_names=TRUE, row_names_gp = gpar(fontsize = 20, fontface="bold"),                           
              cluster_columns=TRUE,                           
              show_column_dend=FALSE,                         
              column_title="Maturation",                      
              column_title_side="top",                          
              column_title_gp=gpar(fontsize=22, fontface="bold"),
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

png("abundant_microbes.png",height=7800,width=15600,res=600)
plot_grid(h1,h2,h3,ncol=3, labels=LETTERS[1:3], rel_widths = c(1, 1)) 
dev.off()
