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