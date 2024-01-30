# Abundant pathway at vegetative stage

v <- read.csv("vegetative_pathways.csv",header=T,row.names=1)
v <- as.matrix(v)
myBreaks <- seq(-3, 3, length.out = 35)
myCol <- colorRampPalette(brewer.pal(3, "Set1"))(35)
heat <- t(scale(t(v)))
tiff("sig_pathway.tiff",height=8800,width=12800,res=600)
Heatmap(heat,                                                 
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
        cluster_columns=TRUE,                           
        show_column_dend=FALSE,                         
        column_title="Pathways abundant at Vegetative Stage",                      
        column_title_side="top",
        column_title_gp=gpar(fontsize=16, fontface="bold"),
        column_title_rot=0,                                
        show_column_names=FALSE,                           
        width = unit(6, "cm"))
dev.off()

# Abundant microbes at the vegetative stage

library(ComplexHeatmap)
library(colorRamp2)
library(RColorBrewer)
data <- read.csv("merged_abundance_table_species1.csv", header=T, row.names=1)  
data = as.matrix(data)  
data2 <- log2(data + 1)
data3 <- read.table("merged_abundance_table_species2.txt",sep=",",header=T,row.names=1)
data3 <- as.matrix(data3)
data3 <- log2(data3 + 1)
data4 <- read.table("merged_abundance_table_species4.txt",sep=",",header=T,row.names=1)
data4 <- as.matrix(data4)
data4 <- log2(data4 + 1)
myBreaks <- seq(-3, 3, length.out = 35)
myCol <- colorRampPalette(brewer.pal(9, "YlOrRd"))(35)
h1 <- Heatmap(data2,
            col=colorRamp2(myBreaks, myCol),
            show_heatmap_legend = FALSE,
              cluster_rows=TRUE,
            show_row_dend=FALSE,
            row_title="Bacterial Species",
            row_title_side="left",
            row_title_gp=gpar(fontsize=16, fontface="bold"),
            row_title_rot=90,
            show_row_names=TRUE,
            row_names_side="left",
            row_names_gp=gpar(fontsize=14),
            cluster_columns=TRUE,
            show_column_dend=FALSE,
            column_title="Medicago truncatula",
            column_title_side="top",
            column_title_gp=gpar(fontsize=16, fontface="bold"),
            column_title_rot=0,
            show_column_names=FALSE, width = unit(6, "cm"))


h2 <- Heatmap(data3, col=colorRamp2(myBreaks, myCol),
            show_heatmap_legend = FALSE,
              cluster_rows=TRUE,
            show_row_dend=FALSE,
            row_title=NA,        
            row_title_side="left",
            row_title_gp=gpar(fontsize=16, fontface="bold"),
            row_title_rot=90,
            show_row_names=TRUE,
            row_names_side="left",
            row_names_gp=gpar(fontsize=14),
            cluster_columns=TRUE,
            show_column_dend=FALSE,
            column_title="Cicer arietinum",
            column_title_side="top",
            column_title_gp=gpar(fontsize=16, fontface="bold"),
            column_title_rot=0,
            show_column_names=FALSE, width = unit(6, "cm"))


h3 <- Heatmap(data4, col=colorRamp2(myBreaks, myCol),
            name="Abundance Z Score", heatmap_legend_param=list(
              color_bar="continuous",
              legend_direction="vertical",
              legend_width=unit(5,"cm"),
              title_position="leftcenter-rot",
              title_gp=gpar(fontsize=15, fontface="bold")),
              cluster_rows=TRUE,
            show_row_dend=FALSE,
            row_title=NA,
            show_row_names=TRUE,
            row_names_side="left",
            row_names_gp=gpar(fontsize=14),
            cluster_columns=TRUE,
            show_column_dend=FALSE,
            column_title="Arachis hypogaea L.",
            column_title_side="top",
            column_title_gp=gpar(fontsize=16, fontface="bold"),
            column_title_rot=0,
            show_column_names=FALSE, width = unit(6, "cm"))
            
library(ggplotify)
library(grid)
 library(cowplot)
library(colorspace)
h1 <- as.ggplot(h1)
 h2 <- as.ggplot(h2)
  h3 <- as.ggplot(h3)
 plot_grid(h1,h2,h3,ncol=3, labels=LETTERS[1:3], rel_widths = c(1, 1)) 
