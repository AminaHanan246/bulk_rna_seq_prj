#==========================================
#Load required Packages
#==========================================
library(tximport)
library(readr)
library(DESeq2)
library(RColorBrewer)
library(gplots)
library(vsn)
library(pheatmap)
library(apeglm)
library(tximportData)
library(limma)
library(sva)
library(tximeta)
library(magrittr)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(vegan)
library(ggplot2)

install.packages("BiocManager")
BiocManager::install(version = "devel")
BiocManager::install("vegan")
getwd()
setwd("D:/SMau/bulkrnaseq_tutorial/fastq")

library(dplyr)
library(tibble)

#==========================================
#set paths and Sample info
#==========================================

#Loading counts after running featureCounts
raw_counts <- read.csv("GSE106305_counts_matrix.csv", header = TRUE, row.names = "Geneid", stringsAsFactors = FALSE)
head(raw_counts)
raw_counts <- raw_counts[,sort(colnames(raw_counts))]
colSums(raw_counts)

#settingnup metadata
condition <- c(rep("LNCAP_Hypoxia", 2), rep("LNCAP_Normoxia", 2), rep("PC3_Hypoxia", 2), rep("PC3_Normoxia", 2))
print(condition)

my_colData <- as.data.frame(condition)
rownames(my_colData) <- colnames(raw_counts)
head(my_colData)

#==========================================
# DESeq - for differential expressed genes
#==========================================

#Create DESeq2 object
dds <- DESeqDataSetFromMatrix(countData = raw_counts,
                              colData = my_colData,
                              design = ~condition)
dds <- DESeq(dds)
dds
counts(dds, normalized = F) #raw counts
normalized_counts <- counts(dds, normalized = T)
head(normalized_counts)

#==========================================
#Ensmeble IDs to gene symbols using BioMArt
#==========================================
