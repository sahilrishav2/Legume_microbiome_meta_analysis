library(microbiome)
r <- read.csv("merged_genera_veg.csv",row.names=1)
transform <- microbiome::transform
rel <- transform(r, "compositional")
meta <- read.csv("metadata.csv")
library(SIAMCAT)
label.dev <- create.label(meta=meta,
                          label='stage', case='Vegetative')
dim(rel)
legume.obj <- siamcat(feat=rel,
                      label=label.dev,
                      meta=meta)
legume.obj.1 <- filter.features(legume.obj,
                                filter.method = 'abundance',
                                cutoff = 0.001)
legume.obj.2 <- check.associations(legume.obj.1, log.n0 = 1e-06, alpha = 0.01)
legume.obj.2 <- normalize.features(legume.obj.2, norm.method = "log.unit",               
                             norm.param = list(log.n0 = 1e-06, n.p = 2,norm.margin = 1))
legume.obj.2 <-  create.data.split(legume.obj.2, num.folds = 5, num.resample = 2) # splitting features for cross validation
legume.obj.2 <- train.model(legume.obj.2, method = "lasso") # training model
legume.obj.2 <- make.predictions(legume.obj.2) # making predictions
pred_matrix <- pred_matrix(legume.obj.2)
legume.obj.2 <-  evaluate.predictions(legume.obj.2) # evaluating predictions
v=legume.obj.2@associations[["assoc.results"]]
write.csv(v,"diff_exp_veg_rep_metagenome_data.csv")

library(ComplexHeatmap)
library(colorRamp2)
library(RColorBrewer)
write.csv(data, "diff_exp_veg_rep_metagenome_data.csv")
dim(data)
v <- as.matrix(data)
myBreaks <- seq(-3, 3, length.out = 35)
myCol <- colorRampPalette(brewer.pal(3, "Set1"))(35)
heat <- t(scale(t(v)))
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
        show_row_dend=TRUE,
        row_title= NA,                               
        show_row_names=TRUE,
        row_names_side="left",                            
        cluster_columns=TRUE,                           
        show_column_dend=FALSE,                         
        #column_title="Pathways abundant at Vegetative Stage",                      
        #column_title_side="top",
        #column_title_gp=gpar(fontsize=16, fontface="bold"),
        #column_title_rot=0,                                
        show_column_names=FALSE,                           
        width = unit(6, "cm"))
