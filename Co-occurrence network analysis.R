########## Co-occurrence network analysis of marker bacterial genera at each developmental stage

#Vegetative

count <- read.csv("countdata.csv",header=T,row.names=1)
meta <- read.csv("metadata.csv",header=T,row.names=1)
dim(count)
dim(meta)
identical(colnames(count),rownames(meta))
rel <- read.csv("veg_batch_adj_p_0.05.csv",header=T,row.names=1)
dim(rel)
library(igraph)   # For network analysis
library(vegan)    # For calculating dissimilarity metrics
library(dplyr)
keep <- rownames(rel[1:175,])
data_subset <- count[rownames(count) %in% keep, ]
dim(data_subset)
developmental_stage <- "Vegetative"  # Specify the developmental stage of interest
meta_stage <- subset(meta, Developmental_stage == developmental_stage)
count_stage <- data_subset[, colnames(data_subset) %in% rownames(meta_stage)]               
dim(count_stage)
write.csv(count_stage,"vegetative_count.csv")
diss_matrix <- vegdist(count_stage, method = "jaccard")
diss_matrix[is.na(diss_matrix)] <- 0
diss_matrix <- as.matrix(diss_matrix)
sim_matrix <- 1 - diss_matrix
sim_vector <- sim_matrix[lower.tri(sim_matrix)]
sorted_sim_vector <- sort(sim_vector, decreasing = TRUE)
sparsity_percentile <- 10
threshold <- quantile(sorted_sim_vector, 1 - sparsity_percentile / 100)  #calculating threshold value
threshold
threshold <- 0.3777005
adj_matrix <- sim_matrix >= threshold
network <- graph.adjacency(adj_matrix, mode = "undirected")
write.graph(network, file = "veg_bacteria_network.graphml", format = "graphml")

#Reproductive

rel <- read.csv("rep_batch_adj_p_0.05.csv",header=T,row.names=1)
dim(rel)
keep <- rownames(rel[1:329,])
data_subset <- count[rownames(count) %in% keep, ]
dim(data_subset)
developmental_stage <- "Reproductive"  # Specify the developmental stage of interest
meta_stage <- subset(meta, Developmental_stage == developmental_stage)
count_stage <- data_subset[, colnames(data_subset) %in% rownames(meta_stage)]               
dim(count_stage)
write.csv(count_stage,"reproductive_count.csv")
diss_matrix <- vegdist(count_stage, method = "jaccard")
diss_matrix[is.na(diss_matrix)] <- 0
diss_matrix <- as.matrix(diss_matrix)
sim_matrix <- 1 - diss_matrix
sim_vector <- sim_matrix[lower.tri(sim_matrix)]
sorted_sim_vector <- sort(sim_vector, decreasing = TRUE)
sparsity_percentile <- 10
threshold <- quantile(sorted_sim_vector, 1 - sparsity_percentile / 100)  #calculating threshold value
threshold
threshold <- 0.4263665
adj_matrix <- sim_matrix >= threshold
network <- graph.adjacency(adj_matrix, mode = "undirected")
write.graph(network, file = "rep_bacteria_network.graphml", format = "graphml")

#Maturation

rel <- read.csv("mat_batch_adj_p_0.05.csv",header=T,row.names=1)
dim(rel)
keep <- rownames(rel[1:143,])

data_subset <- count[rownames(count) %in% keep, ]
dim(data_subset)
developmental_stage <- "Maturation"  # Specify the developmental stage of interest                                                                
meta_stage <- subset(meta, Developmental_stage == developmental_stage)
count_stage <- data_subset[, colnames(data_subset) %in% rownames(meta_stage)] 
dim(count_stage)

write.csv(count_stage,"maturation_count.csv")
diss_matrix <- vegdist(count_stage, method = "jaccard")
diss_matrix[is.na(diss_matrix)] <- 0
diss_matrix <- as.matrix(diss_matrix)
sim_matrix <- 1 - diss_matrix
sim_vector <- sim_matrix[lower.tri(sim_matrix)]
sorted_sim_vector <- sort(sim_vector, decreasing = TRUE)
sparsity_percentile <- 10
threshold <- quantile(sorted_sim_vector, 1 - sparsity_percentile / 100)  #calculating threshold value
threshold
threshold <- 0.391506
adj_matrix <- sim_matrix >= threshold
network <- graph.adjacency(adj_matrix, mode = "undirected")
write.graph(network, file = "mat_bacteria_network.graphml", format = "graphml")