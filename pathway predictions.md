## For functional predictions using picrust2
conda activate picrust2
picrust2_pipeline.py -s asv_seqs.fna -i asv_counts.tsv -o picrust2_out_pipeline -p 15

## Loading R and performing differential abundance analysis of merged pathways from all studies

# For vegetative stage

library(edgeR)
p <- read.csv("pathway.csv",header=T,row.names=1)
dim(p)
d0 <- DGEList(p)
d0 <- calcNormFactors(d0)
cutoff <- 1
drop <- which(apply(cpm(d0), 1, max) < cutoff)
d <- d0[-drop,] 
dim(d)
m = read.csv("meta_data.csv", header=T)
identical(colnames(p),m$Sample)
identical(colnames(d),m$Sample)
dim(m)
head(m)
snames <-  m$Developmental_stage
head(snames)
group <- snames
mm <- model.matrix(~0 + group)
y <- voom(d, mm, plot = T)
tmp <- voom(d0, mm, plot = T)
fit <- lmFit(y, mm)
head(coef(fit))
contr <- makeContrasts(groupVegetative - groupOthers, levels = colnames(coef(fit)))
tmp <- contrasts.fit(fit, contr)
tmp <- eBayes(tmp)
top.table <- topTable(tmp, sort.by = "P", n = Inf)
head(top.table, 20)
length(which(top.table$adj.P.Val < 0.05))
newdf <- top.table[top.table$adj.P.Val < 0.05,]
dim(newdf)
newdf$Pathway <- rownames(newdf)
newdf <- newdf[,c("Pathway", names(newdf)[1:6])]
write.csv(newdf,"diff_veg.csv")

# For reproductive stage

m2 = read.csv("meta_data2.csv", header=T)
identical(colnames(p),m$Sample)
identical(colnames(d),m$Sample)
dim(m)
head(m)
snames <-  m$Developmental_stage
head(snames)
group <- snames
mm <- model.matrix(~0 + group)
y <- voom(d, mm, plot = T)
tmp <- voom(d0, mm, plot = T)
fit <- lmFit(y, mm)
head(coef(fit))
contr <- makeContrasts(groupReproductive - groupOthers, levels = colnames(coef(fit)))
tmp <- contrasts.fit(fit, contr)
tmp <- eBayes(tmp)
top.table <- topTable(tmp, sort.by = "P", n = Inf)
head(top.table, 20)
length(which(top.table$adj.P.Val < 0.05))
newdf <- top.table[top.table$adj.P.Val < 0.05,]
dim(newdf)
newdf$Pathway <- rownames(newdf)
newdf <- newdf[,c("Pathway", names(newdf)[1:6])]
write.csv(newdf,"diff_rep.csv")

# For maturation stage

m3 = read.csv("meta_data3.csv", header=T)
identical(colnames(p),m$Sample)
identical(colnames(d),m$Sample)
dim(m)
head(m)
snames <-  m$Developmental_stage
head(snames)
group <- snames
mm <- model.matrix(~0 + group)
y <- voom(d, mm, plot = T)
tmp <- voom(d0, mm, plot = T)
fit <- lmFit(y, mm)
head(coef(fit))
contr <- makeContrasts(groupMaturation - groupOthers, levels = colnames(coef(fit)))
tmp <- contrasts.fit(fit, contr)
tmp <- eBayes(tmp)
top.table <- topTable(tmp, sort.by = "P", n = Inf)
head(top.table, 20)
length(which(top.table$adj.P.Val < 0.05))
newdf <- top.table[top.table$adj.P.Val < 0.05,]
dim(newdf)
newdf$Pathway <- rownames(newdf)
newdf <- newdf[,c("Pathway", names(newdf)[1:6])]
write.csv(newdf,"diff_mat.csv")