## To perform ASV based analysis, we directly load paired-end reads fastq files in R and perform processing through dada2

library(dada2)
path <- "/home/sahil/Documents/Legume/Amplicon1" #directory containing fastq files
list.files(path)

## Sorting of fastq files

fnFs <- sort(list.files(path, pattern="_1.fastq", full.names=TRUE))
fnRs <- sort(list.files(path, pattern="_2.fastq", full.names=TRUE))

## Extraction of sample names having extension ".fastq"

sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

## Inspecting read quality profiles

plotQualityProfile(fnFs[1:2]) #forward reads
plotQualityProfile(fnRs[1:2]) #reverse reads

## Assigning the filenames for the filtered fastq.gz files

filtFs <- paste0(sample.names, "-F-filt.fastq.gz")
filtRs <- paste0(sample.names, "-R-filt.fastq.gz")

## For filtering of reads

for(i in seq_along(fnFs)) { fastqPairedFilter(c(fnFs[i], fnRs[i]), c(filtFs[i], filtRs[i]), truncLen=c("#forward reads truncate position","#reverse reads truncate position"), maxN=0, maxEE=c(2,2), truncQ=2, compress=TRUE, verbose=TRUE) }

## Sample inference

derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)
names(derepFs) <- sample.names
names(derepRs) <- sample.names
dadaFs <- dada(derepFs, err=inflateErr(tperr1,3), selfConsist = TRUE, MAX_CONSIST=20)
dadaRs <- dada(derepRs, err=inflateErr(tperr1,3), selfConsist = TRUE, MAX_CONSIST=20)

## Merging of paired reads

mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)
head(mergers[[1]])

## Constructing sequence table

seqtab <- makeSequenceTable(mergers)
dim(seqtab)
table(nchar(getSequences(seqtab))) #Inspecting distribution of sequence lengths

## Removing chimeras

seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab)

## Assigning taxonomy

taxa <- assignTaxonomy(seqtab.nochim, "silva_nr99_v138.1_wSpecies_train_set.fa.gz",multithread=TRUE)

## Inspecting the taxonomic assignments

taxa.print <- taxa
rownames(taxa.print) <- NULL
head(taxa.print)

## Extracting the count values and the sequences from dada2 output

# giving our seq headers more manageable names (ASV_1, ASV_2...)
asv_seqs <- colnames(seqtab.nochim)
asv_headers <- vector(dim(seqtab.nochim)[2], mode="character")

for (i in 1:dim(seqtab.nochim)[2]) {
  asv_headers[i] <- paste(">ASV", i, sep="_")
}

# making and writing out a fasta of our final ASV seqs:
asv_fasta <- c(rbind(asv_headers, asv_seqs))
write(asv_fasta, "ASVs.fa")

# count table:
asv_tab <- t(seqtab.nochim)
row.names(asv_tab) <- sub(">", "", asv_headers)
write.table(asv_tab, "ASVs_counts.tsv", sep="\t", quote=F, col.names=NA)

# taxonomy table:
rownames(taxa.print) <- rownames(asv_tab)
head(taxa.print)
write.table(taxa.print, "ASVs_taxonomy.tsv", sep="\t", quote=F, col.names=NA)

## Now, we load phyloseq for further analysis of the amplicon data

library(phyloseq)
seq <- as.matrix(seqtab.nochim)+1 # for removing zeroes from the count table
tax <- taxa
dim(seq)
colnames(seq) <- paste0("ASV_", "#no. of ASV ids") #for giving unique ASV Ids to each sequence
rownames(tax) <- colnames(seq)
ps <- phyloseq(otu_table(seq, taxa_are_rows=FALSE), 
               tax_table(tax))
library(microbiomeutilities)
ps1 <- ps #assigning all the information of ps to a new object ps1
rank_names(ps1) # check the taxonomy levels
#Removing sequences not aligning to the bacteria
ps1a <- subset_taxa(ps1,Order!="Chloroplast")
ps1b <- subset_taxa(ps1a,Family!="Mitochondria")
ps1 <- subset_taxa(ps1b,Kingdom!="Archaea")
ps
ps1
write.table(ps1@otu_table, "ps1_ASVs_counts.tsv", sep="\t", quote=F, col.names=NA) ## writing a new count file
write.table(ps1@tax_table, "ps1_ASVs_taxonomy.tsv", sep="\t", quote=F, col.names=NA) ## writing a new taxonomy file
