############ Corrrelation between marker bacteria and abundant microbial pathways of each developmental stage
## spearmann correlation test:
  
library(corrplot)
library(psych)
p <- read.csv("path.csv",row.names=1)
s <- read.csv("species.csv",row.names=1)
m1 <-corr.test(p,s,method="spearman",adjust="fdr")
m2 <- data.frame(m1["p.adj"])
matrix_corr <- cor(p,s,method="spearman",use="complete.obs")
colnames(m2) <- colnames(matrix_rel_eth)
m2 <- as.matrix(m2)
str(matrix_corr)
str(m2)
png("corrplot.png",height=5400,width=9800,res=600)
corrplot(matrix_corr,
         method = "circle",
         outline = T, 
         addgrid.col = "grey", 
         p.mat= matrix_q,
         sig.level = 0.01,
         insig ="blank", 
         addrect = 4, 
         cl.pos = "b", 
         tl.col = "black", #colour of text labels 
         tl.cex = 0.5, #Size of text labels
         cl.cex = 0.5)
dev.off()
