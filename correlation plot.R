############ Corrrelation between marker bacteria and abundant microbial pathways of each developmental stage
## spearmann correlation test:
  
library(corrplot)
library(psych)
p <- read.csv("path3.csv",row.names=1)
s <- read.csv("species3.csv",row.names=1)
metrix_p <-corr.test(p,s,method="spearman",adjust="fdr")
matrix_q <- data.frame(metrix_p["p.adj"])
matrix_rel_eth <- cor(p,s,method="spearman",use="complete.obs")
colnames(matrix_q) <- colnames(matrix_rel_eth)
matrix_q <- as.matrix(matrix_q)
str(matrix_rel_eth)
str(matrix_q)
png("corr_mat2.png",height=5400,width=9800,res=600)
corrplot(matrix_rel_eth,
         #title = "Correlation Plot",
         method = "circle",
         outline = T, 
         addgrid.col = "grey", 
         p.mat= matrix_q,
         sig.level = 0.01,
         insig ="blank",
         #mar = c(0,0,1,0), 
         addrect = 4, 
         #rect.col = "grey", 
         #rect.lwd = 1,
         cl.pos = "b", 
         #        tl.pos = "ld",
         tl.col = "black", #colour of text labels 
         tl.cex = 0.5, #Size of text labels
         cl.cex = 0.5)
dev.off()