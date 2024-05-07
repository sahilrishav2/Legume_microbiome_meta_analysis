###### Bacterial community analyses
## Alpha diversity analysis of studies

library(phyloseq)
library(ggplot2) 
library(RColorBrewer)
count <- read.csv("countdata.csv",header=T,row.names=1)
meta <- read.csv("metadata.csv",header=T,row.names=1)
ps <- phyloseq(otu_table(count, taxa_are_rows=TRUE), 
               sample_data(meta))
ps
p <- plot_richness(ps, "Study", color="Developmental_stage", measures = c("Shannon","Simpson"))
p + geom_boxplot(aes(x=Study, colour=Developmental_stage)) +
  theme_bw() +                             
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        panel.background = element_blank(), panel.border = element_rect(colour = "black")) + theme(axis.text.x = element_text(hjust = 1, size = 16, face = "bold"),  axis.title.x = element_text(size = 20), # Increase x-axis label font size
                                                                                                   axis.title.y = element_text(size = 18), # Increase y-axis label font size
                                                                                                   legend.title=element_text(size=18), 
                                                                                                   legend.text=element_text(size=16)) +  coord_flip() +   theme(axis.text.y = element_text(size = 16, face = "bold"))

## Between class PCA for samples of different legumes

# For plotting pca plot:

library(ade4)
library("factoextra")
library(dplyr)
rel <- read.csv("countdata.csv",header=T,row.names=1)
rel2 <- as.matrix(t(rel))
meta <- read.csv("metadata.csv",header=T)
res.pca <- dudi.pca(rel2,
                    scannf = FALSE, nf = 5)
groups <- as.factor(meta$Plant)
s.class(res.pca$li,                                                  
        fac = groups,  # color by groups
        col = c("brown4",  "lightsalmon2", "darkslategray1", "blue3", "gold2", "firebrick1", "lightseagreen", "darkmagenta"), clabel = 0  
)


rand1 <- randtest(bca(res.pca, groups, scan = FALSE), 999)
rand1
