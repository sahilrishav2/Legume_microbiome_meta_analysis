##### Microbial composition of included datasets

library(ggplot2)
library(reshape2)
df <- read.table("phylum_study_abundance_2.csv", header = TRUE, row.names = 1, sep="\t")
dim(df)
df$Phylum <- rownames(df)
rownames(df) <- NULL
melted_df <- melt(df, id.vars = "Phylum")
melted_df$Phylum <- factor(melted_df$Phylum, levels = df$Phylum)
plot <- ggplot(melted_df, aes(x = variable, y = value, fill = Phylum)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(x = "Study", y = "Mean Relative Abundance", title = "Phylum Abundance in Different Studies") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 65, hjust = 1, size = 18, face = "bold"),  axis.title.x = element_text(size = 20, face = "bold"), # Increase x-axis label font size
        axis.title.y = element_text(size = 20, face = "bold"), # Increase y-axis label font size
        plot.title = element_text(size = 24, face = "bold"), legend.title=element_text(size=18, face = "bold"), 
        legend.text=element_text(size=16, face = "bold")) +  theme(axis.text.y = element_text(size = 16, face = "bold")) +
  scale_fill_brewer(palette = "Set3")


plot


##### Alluvial plot showing the association of studies and developmental stages with the microbial community at the phylum level

library(ggalluvial)
library(reshape2)
b <- read.csv("phylum_sample_wise.csv",header=T)
dim(b)
abundance_melt <- melt(b, id.vars="Phylum")
m2 <- read.csv("meta4.csv",header=T)
data2 <- merge(abundance_melt, m2, by.x="variable", by.y="Sample")
library(RColorBrewer)
png("alluvial_plot.png",height=7400,width=10800,res=600)
ggplot(data = data2,
       aes(axis1 = Study, axis2 = Developmental_stage, axis3 = Phylum, y = value)) +
  geom_alluvium(aes(fill = Phylum)) +
  geom_stratum() + geom_text(stat = "stratum",
                             aes(label = after_stat(stratum))) +
  scale_x_discrete(limits = c("Developmental_stage", "Phylum", "Study"),
                   expand = c(0.25, 0.1)) + scale_fill_brewer(palette = "Paired") +
  theme_void() 
dev.off()