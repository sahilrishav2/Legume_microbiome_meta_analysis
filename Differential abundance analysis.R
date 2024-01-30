##### Batch effect adjustment
m <- read.csv("metadata.csv",header=T,row.names=1)
r <- read.csv("abundancedata.csv",header=T,row.names=1)
identical(rownames(m),colnames(r))
library(microbiome)
transform <- microbiome::transform
rel <- transform(r, "compositional") ##### converting abundance values into relative abundances
dim(rel)
library(MMUPHin)
meta <- m
dim(m)
meta$Study <- as.factor(meta$Study)
fit_adjust_batch <- adjust_batch(feature_abd = rel,
                                 batch = "Study",
                                 covariates = "Developmental_stage",
                                 data = meta,
                                 control = list(verbose = FALSE))
f <- fit_adjust_batch$feature_abd_adj
write.csv(f,"countdata.csv")

###### Differential abundance analysis
## Vegetative stage
rel <- read.csv("countdata.csv",header=T,row.names=1)
library(SIAMCAT)
label.dev <- create.label(meta=meta,
                          label='Developmental_stage', case='Vegetative')

legume.obj <- siamcat(feat=rel,
                      label=label.dev,
                      meta=meta)
legume.obj.1 <- filter.features(legume.obj,
                                filter.method = 'abundance',
                                cutoff = 0.001)
legume.obj.2 <- check.associations(legume.obj.1, log.n0 = 1e-06, alpha = 0.01)
legume.obj.2 <- normalize.features(legume.obj.2, norm.method = "log.unit",               
                             norm.param = list(log.n0 = 1e-06, n.p = 2,norm.margin = 1)) # normalizing features
legume.obj.2 <-  create.data.split(legume.obj.2, num.folds = 5, num.resample = 2) # splitting features for cross validation
legume.obj.2 <- train.model(legume.obj.2, method = "lasso") # training model
legume.obj.2 <- make.predictions(legume.obj.2) # making predictions
pred_matrix <- pred_matrix(legume.obj.2)
legume.obj.2 <-  evaluate.predictions(legume.obj.2) # evaluating predictions
v=legume.obj.2@associations[["assoc.results"]]
write.csv(v,"batch_adj_v_asso.csv")

## Reproductive stage
label.dev <- create.label(meta=meta,
                          label='Developmental_stage', case='Reproductive')

legume.obj.3 <- siamcat(feat=rel,
                        label=label.dev,
                        meta=meta)
legume.obj.4 <- filter.features(legume.obj.3,
                                filter.method = 'abundance',
                                cutoff = 0.001)
legume.obj.5 <- check.associations(legume.obj.4, log.n0 = 1e-06, alpha = 0.01)
legume.obj.5 <- normalize.features(legume.obj.5, norm.method = "log.unit",               
                                   norm.param = list(log.n0 = 1e-06, n.p = 2,norm.margin = 1)) # normalizing features
legume.obj.5 <-  create.data.split(legume.obj.5, num.folds = 5, num.resample = 2) # splitting features for cross validation
legume.obj.5 <- train.model(legume.obj.5, method = "lasso") # training model
legume.obj.5 <- make.predictions(legume.obj.5) # making predictions
pred_matrix <- pred_matrix(legume.obj.5)
legume.obj.5 <-  evaluate.predictions(legume.obj.5) # evaluating predictions
r <- legume.obj.5@associations[["assoc.results"]]
write.csv(r,"batch_adj_r_asso.csv")

## Maturation stage
label.dev <- create.label(meta=meta,
                          label='Developmental_stage', case='Maturation')

legume.obj.6 <- siamcat(feat=rel,
                        label=label.dev,
                        meta=meta)
legume.obj.7 <- filter.features(legume.obj.6,
                                filter.method = 'abundance',
                                cutoff = 0.001)
legume.obj.8 <- check.associations(legume.obj.7, log.n0 = 1e-06, alpha = 0.01)
legume.obj.8 <- normalize.features(legume.obj.8, norm.method = "log.unit",               
                                   norm.param = list(log.n0 = 1e-06, n.p = 2,norm.margin = 1)) # normalizing features
legume.obj.8 <-  create.data.split(legume.obj.8, num.folds = 5, num.resample = 2) # splitting features for cross validation
legume.obj.8 <- train.model(legume.obj.8, method = "lasso") # training model
legume.obj.8 <- make.predictions(legume.obj.8) # making predictions
pred_matrix <- pred_matrix(legume.obj.8)
legume.obj.8 <-  evaluate.predictions(legume.obj.8) # evaluating predictions
m=legume.obj.8@associations[["assoc.results"]]
write.csv(m,"batch_adj_m_asso.csv")