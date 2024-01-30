## Enriched phyla at each developmental stage

library(ggplot2)
library(RColorBrewer)

## Reproductive stage 
r <- read.csv("rep_batch_adj_p_0.05.csv", header = TRUE)

# Create a data frame with counts of each Phylum
phylum_counts <- as.data.frame(table(r$Phylum))
phylum_counts <- phylum_counts[order(phylum_counts$Freq, decreasing = TRUE), ]

manual_colors <- c(
  "Proteobacteria" = "red",
  "Actinobacteriota" = "blue",
  "Firmicutes" = "green",
  "Bacteroidota" = "purple",
  "Verrucomicrobiota" = "orange",
  "Planctomycetota" = "pink",
  "Myxococcota" = "brown",
  "Acidobacteriota" = "gray",
  "Chloroflexi" = "cyan",
  "Cyanobacteria" = "yellow",
  "Desulfobacterota" = "magenta",
  "Armatimonadota" = "darkgreen",
  "Bdellovibrionota" = "violet",
  "Deinococcota" = "darkred",
  "Gemmatimonadota" = "darkblue",
  "Spirochaetota" = "lightblue",
  "Abditibacteriota" = "lightgreen",
  "Elusimicrobiota" = "lightgray",
  "Entotheonellaeota" = "darkorange",
  "Fibrobacterota" = "paleturquoise1",
  "Fusobacteriota" = "orchid2",
  "Nitrospirota" = "gold",
  "Sumerlaeota" = "lightsalmon2"
)

data <- data.frame(
  Phylum = factor(r$Phylum, levels = phylum_counts$Var1),  # Reorder Phylum based on abundance
  Adjusted_Pvalue = r$p.adj,
  AUROC = r$auc
)

# Filter the data based on the AUROC cutoff
filtered_data <- subset(data, AUROC >= 0.5)

# Create the plot
plot <- ggplot(data, aes(x = AUROC, y = -log10(Adjusted_Pvalue), color = Phylum)) +
  geom_point(size = 3, color = "black") +
  geom_vline(xintercept = 0.5, linetype = "dashed", color = "black", size = 0.5) +
  geom_point(data = filtered_data, size = 6) +
  geom_text(aes(x = 0.2, y = 10, label = "Depleted in Reproductive stage"), color = "black", size = 4, hjust = 0) +
  geom_text(aes(x = 0.8, y = 10, label = "Enriched in Reproductive stage"), color = "red", size = 4, hjust = 1) +
  scale_color_manual(values = manual_colors) +  # Apply manual colors
  labs(
    x = "AUROC",
    y = "-log10(Adjusted P-value)",
    color = "Phylum"
  ) +
  theme_minimal()

# Display the plot
print(plot)


## Vegetative stage
v <- read.csv("veg_batch_adj_p_0.05.csv", header = TRUE)

# Create a data frame with counts of each Phylum
phylum_counts <- as.data.frame(table(v$Phylum))
phylum_counts <- phylum_counts[order(phylum_counts$Freq, decreasing = TRUE), ]

manual_colors <- c(
  "Proteobacteria" = "red",
  "Actinobacteriota" = "blue",
  "Firmicutes" = "green",
  "Bacteroidota" = "purple",
  "Verrucomicrobiota" = "orange",
  "Planctomycetota" = "pink",
  "Myxococcota" = "brown",
  "Acidobacteriota" = "gray",
  "Chloroflexi" = "cyan",
  "Cyanobacteria" = "yellow",
  "Desulfobacterota" = "magenta",
  "Armatimonadota" = "darkgreen",
  "Bdellovibrionota" = "violet",
  "Deinococcota" = "darkred",
  "Gemmatimonadota" = "darkblue",
  "Spirochaetota" = "lightblue",
  "Abditibacteriota" = "lightgreen",
  "Elusimicrobiota" = "lightgray",
  "Entotheonellaeota" = "darkorange",
  "Fibrobacterota" = "paleturquoise1",
  "Fusobacteriota" = "orchid2",
  "Nitrospirota" = "gold",
  "Sumerlaeota" = "lightsalmon2"
)

data <- data.frame(
  Phylum = factor(v$Phylum, levels = phylum_counts$Var1),  # Reorder Phylum based on abundance
  Adjusted_Pvalue = v$p.adj,
  AUROC = v$auc
)

# Filter the data based on the AUROC cutoff
filtered_data <- subset(data, AUROC >= 0.5)

# Create the plot
plot <- ggplot(data, aes(x = AUROC, y = -log10(Adjusted_Pvalue), color = Phylum)) +
  geom_point(size = 3, color = "black") +
  geom_vline(xintercept = 0.5, linetype = "dashed", color = "black", size = 0.5) +
  geom_point(data = filtered_data, size = 6) +
  geom_text(aes(x = 0.2, y = 10, label = "Depleted in Vegetative stage"), color = "black", size = 4, hjust = 0) +
  geom_text(aes(x = 0.8, y = 10, label = "Enriched in Vegetative stage"), color = "red", size = 4, hjust = 1) +
  scale_color_manual(values = manual_colors) +  # Apply manual colors
  labs(
    x = "AUROC",
    y = "-log10(Adjusted P-value)",
    color = "Phylum"
  ) +
  theme_minimal()

# Display the plot
print(plot)



## Maturation stage
m <- read.csv("mat_batch_adj_p_0.05.csv", header = TRUE)

# Create a data frame with counts of each Phylum
phylum_counts <- as.data.frame(table(m$Phylum))
phylum_counts <- phylum_counts[order(phylum_counts$Freq, decreasing = TRUE), ]

manual_colors <- c(
  "Proteobacteria" = "red",
  "Actinobacteriota" = "blue",
  "Firmicutes" = "green",
  "Bacteroidota" = "purple",
  "Verrucomicrobiota" = "orange",
  "Planctomycetota" = "pink",
  "Myxococcota" = "brown",
  "Acidobacteriota" = "gray",
  "Chloroflexi" = "cyan",
  "Cyanobacteria" = "yellow",
  "Desulfobacterota" = "magenta",
  "Armatimonadota" = "darkgreen",
  "Bdellovibrionota" = "violet",
  "Deinococcota" = "darkred",
  "Gemmatimonadota" = "darkblue",
  "Spirochaetota" = "lightblue",
  "Abditibacteriota" = "lightgreen",
  "Elusimicrobiota" = "lightgray",
  "Entotheonellaeota" = "darkorange",
  "Fibrobacterota" = "paleturquoise1",
  "Fusobacteriota" = "orchid2",
  "Nitrospirota" = "gold",
  "Sumerlaeota" = "lightsalmon2"
)

data <- data.frame(
  Phylum = factor(m$Phylum, levels = phylum_counts$Var1),  # Reorder Phylum based on abundance
  Adjusted_Pvalue = m$p.adj,
  AUROC = m$auc
)

# Filter the data based on the AUROC cutoff
filtered_data <- subset(data, AUROC >= 0.5)

# Create the plot
plot <- ggplot(data, aes(x = AUROC, y = -log10(Adjusted_Pvalue), color = Phylum)) +
  geom_point(size = 3, color = "black") +
  geom_vline(xintercept = 0.5, linetype = "dashed", color = "black", size = 0.5) +
  geom_point(data = filtered_data, size = 6) +
  geom_text(aes(x = 0.2, y = 10, label = "Depleted in maturation stage"), color = "black", size = 4, hjust = 0) +
  geom_text(aes(x = 0.8, y = 10, label = "Enriched in maturation stage"), color = "red", size = 4, hjust = 1) +
  scale_color_manual(values = manual_colors) +  # Apply manual colors
  labs(
    x = "AUROC",
    y = "-log10(Adjusted P-value)",
    color = "Phylum"
  ) +
  theme_minimal()

# Display the plot
print(plot)