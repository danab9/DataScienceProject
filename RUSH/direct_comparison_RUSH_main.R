library(DESeq2)
library(readr)
library(stringr)
library(dplyr)
library(ggplot2)


# Import count tables
##---------------------------------------------------------------------------------------##

# Rush Data

# female AD case, 29 cases and 39672 genes
rush_AD <- read.csv("~/DataScienceProject/RUSH/data/all_counts_female.csv",  header=TRUE)
rush_AD <- rush_AD[,2:ncol(rush_AD)]
rush_AD<- dplyr::select(rush_AD, -c(ensembl_gene_id, gene_biotype))
rownames(rush_AD)<- rush_AD$hgnc_symbol
rush_AD <- dplyr::select(rush_AD, -hgnc_symbol)



# female control cases, 13 samples and 39670 genes
rush_control <- read.csv("~/DataScienceProject/RUSH/data/all_counts_control.csv",  header=TRUE)
rush_control <- rush_control[,2:ncol(rush_control)]
rush_control <- dplyr::select(rush_control, -c(ensembl_gene_id, gene_biotype))
rownames(rush_control)<- rush_control$hgnc_symbol
rush_control <- dplyr::select(rush_control, -hgnc_symbol)


# AD and control samples from main paper, 27135 different genes and 30 patients

maindata <- read.table("~/DataScienceProject/RNA-seq/GSE153873_summary_count.star.txt", sep = "\t",header = TRUE,row.names = 1)
colnames(maindata) <- paste0(data.frame(str_split(colnames(maindata),'[.]'))[3,],data.frame(str_split(colnames(maindata),'[.]'))[1,])
dim(maindata)

# Groups
AD <- dplyr::select(maindata, starts_with("AD"))
young <- dplyr::select(maindata, starts_with("Young"))
old <- dplyr::select(maindata, starts_with("Old"))


# Differential gene expression: 
##---------------------------------------------------------------------------------------##

# merge data: 
# control (main old  + rush ) vs AD (main + rush)

# for balanced group, reducing sample size of rush 
rush_AD <- rush_AD[, 1:12]
rush_control <- rush_control[,1:10]

# intersect the measured genes from the data sets

total_control<- merge(old, rush_control,  by= 0 )
rownames(total_control)<- total_control$Row.names
total_control<- dplyr::select(total_control, -Row.names)

total_AD<- merge(AD, rush_AD,  by= 0 )
rownames(total_AD)<- total_AD$Row.names
total_AD<- dplyr::select(total_AD, -Row.names)

total_counts<- merge(total_control, total_AD, by=0)
rownames(total_counts)<- total_counts$Row.names
total_counts<- dplyr::select(total_counts, -Row.names)


total_rush<- merge(rush_control, rush_AD, by=0)
rownames(total_rush)<- total_rush$Row.names
total_rush<- dplyr::select(total_rush, -Row.names)


total_main<- merge(old, AD, by=0)
rownames(total_main)<- total_main$Row.names
total_main<- dplyr::select(total_main, -Row.names)

# DESeq2  attempt
batch<- factor(c(rep("main", ncol(old)),rep("rush",ncol(rush_control)),rep("main", ncol(AD)),rep("rush", ncol(rush_AD))))

condition <- factor(c(rep("C", ncol(old)+ ncol(rush_control)),rep("AD", ncol(AD)+ ncol(rush_AD))))

# DeSeq2 modeling with batch effect
dds <- DESeqDataSetFromMatrix(countData = total_counts,
                              DataFrame(condition,batch), ~ batch+condition)

dds.result <- DESeq(dds)
dds.result <- results(dds.result ,alpha=0.05)
summary(dds.result)

pvalues <- p.adjust(dds.result$pvalue, method = "BH")
sum(pvalues<0.05,na.rm = T)

#-----------------------------------------------------------------------------------------#
# Data Exploration 
## Investigation of Batch Effects 
# main: main samples, rush: rush samples

dds$batch <- batch

# Principal Component Analysis

vsd<- vst(dds, blind=FALSE)

# plot PCA without batch removal 

pcaData <- plotPCA(vsd, intgroup=c("condition", "batch"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=condition, shape = batch)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) 

#ggsave("pca_before_batch_removal.png")

# PCA with batch removal 
assay(vsd) <- limma::removeBatchEffect(assay(vsd), vsd$batch)
pcaData <- plotPCA(vsd, intgroup=c("condition", "batch"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=condition, shape = batch)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) 

#ggsave("PCA_after_batch_removal.png")


# using combat-seq to adjust rush RNA-seq data

library(sva)
library(edgeR)

combatseq_counts<- ComBat_seq(as.matrix(total_counts), batch = batch, group=condition )

# PCA values
pca_combatseq<- as.data.frame(prcomp(combatseq_counts)[2]$rotation)
pca_combatseq$condition<- condition 
pca_combatseq$batch<- batch
ggplot(pca_combatseq, aes(PC1, PC2, color=condition, shape = batch)) +
  geom_point(size=3) 

#ggsave("PCA_after_combatseq_batch_correction.png")


# heatmap of samples

library("pheatmap")
library("RColorBrewer")
sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd$condition, vsd$batch, sep="-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pdf("batchcorrection_heatmap.pdf", height = 4, width = 5)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)
dev.off()


### difference between rush and main is just too huge , thus a direct comparison makes no sense
## after batch removal, there is now no difference visible between the two conditions


