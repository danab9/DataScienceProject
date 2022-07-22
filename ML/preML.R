
library(dplyr)

# all genes (normalized)
library(ggplot2)
library(gplots)
library(lattice)

# load chipseq data
# chipseq_12k <- read.csv("ChIP-seq/iCluster_mat_12k.txt", sep='\t',header = FALSE)
# chipseq_12k <- data.frame(t(chipseq_12k))

chipseq_12k <- read.csv("../ChIP-seq/iCluster_mat_input.txt", sep='\t', header=FALSE)  # regions names
chipseq_12k <- data.frame(t(chipseq_12k))
colnames(chipseq_12k) <- chipseq_12k[34,]  # gene name as column name
chipseq_12k <- chipseq_12k[4:30,] # keep only counts data 

rownames(chipseq_12k) <- c("YoungX1", "YoungX3", "YoungX4", "YoungX5", "YoungX6",
                           "YoungX7", "YoungX8", "YoungX9", "OldX10", "OldX11",
                           "OldX12", "OldX13", "OldX14", "OldX15", "OldX16", 
                           "OldX17", "OldX18", "OldX19",
                           "ADX20", "ADX21", "ADX22", "ADX23","ADX25",
                           "ADX26", "ADX27", "ADX28", "ADX31") 


normalized <- data.frame(t(read.csv("results/normalizedRNA_new.tsv")))
colnames(normalized) <- normalized[1,]
normalized <- normalized[-1,]


rnaseq_all <- normalized
rnaseq_all <- rnaseq_all[which(rownames(rnaseq_all) %in% row.names(chipseq_12k)),]  # remove samples that are not in chipseq

rnaseq_all <- rnaseq_all[order(row.names(rnaseq_all)),] # sort same as chipseq
rnaseq_all <- data.matrix(rnaseq_all)
chipseq_12k <- chipseq_12k[which(row.names(chipseq_12k) %in% row.names(rnaseq_all)),]
chipseq_12k <- chipseq_12k[order(row.names(chipseq_12k)),]

