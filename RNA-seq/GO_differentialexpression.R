library(stringr)
library(DESeq2)
library(dplyr)
library(apeglm)
library(readxl)
library(pheatmap)
library(mixOmics)

source('differentialexpression_function.R')

# import data
rnaseqdata <- read.table("RNAseq/GSE153873_summary_count.star.txt", sep = "\t",header = TRUE,row.names = 1)
colnames(rnaseqdata) <- paste0(data.frame(str_split(colnames(rnaseqdata),'[.]'))[3,],data.frame(str_split(colnames(rnaseqdata),'[.]'))[1,])


# normalized rna
rnaseqdata <- read.table("RNAseq/normalizedRNA.txt")

# import genes from GO enrichment "Regulation of transcription"
genes_regtrans <- read.table("RNAseq/GOenrichment/reg_transcription_own.txt", sep = "\t",header = TRUE)
colnames(genes_regtrans) <- c("Gene Name","Species")

genes_regtrans_paper <- read.table("RNAseq/GOenrichment/reg_transcription_paper.txt", sep = "\t",header = TRUE)
colnames(genes_regtrans_paper) <- c("Gene Name","Species")

#import upregulated genes
# tree main genes from the paper
imp <- c( 'EP300', 'CREBBP','TRRAP')
which(rownames(genes_regtrans)%in% imp)

length(intersect(rownames(genes_regtrans),rownames(genes_regtrans_paper)))

# Use only genes from genes_regtrans
#rnaseqdataGenes <- rnaseqdata$refGenerefGene[unlist(rnaseqdata$refGenerefGene)%in%rownames(genes_regtrans_paper)]
rnaseqdata <- rnaseqdata[rownames(rnaseqdata)%in%rownames(genes_regtrans),]
AD <-  dplyr::select(rnaseqdata, starts_with("AD"))
young <-  dplyr::select(rnaseqdata, starts_with("Young"))
old <-  dplyr::select(rnaseqdata, starts_with("Old"))

rnaseqdata <- data.frame(young,old,AD)


# Comparing gene expression between young AD and old samples

dds <- diffexpression(young,old,AD,tidy=TRUE,result=FALSE)
dds

# plotdf<- data_frame(log2FoldChange=dds$log2FoldChange, padj= dds$padj,pval=dds$pvalue, genes=rownames(rnaseqdata))
# plotdf<- plotdf %>%
#   mutate(gene_type = case_when(log2FoldChange > 0 & pval < 0.05 ~ "upregulated",
#                                log2FoldChange < 0 & pval < 0.05 ~ "downregulated",
#                                TRUE ~ "not significant"))   
# print(table(plotdf$gene_type))
# plotdf <- plotdf[order(plotdf$log2FoldChange,decreasing = TRUE),]
#Heatmap of the count matrix + cluster

ntd <- assay(normTransform(dds))

library("gplots")
png(file = "RT_heatmap.png", height = 1800) 
heatmap.2(ntd, scale = "none", col = bluered(200),keysize=0.75, key.par = list(cex=0.3),key=FALSE,
          trace = "none", density.info = "none",dendrogram = "row",Colv = FALSE, cexRow=0.9)
dev.off()

library("gplots")
png(file = "RT_heatmapNew.png", height = 1800)
heatmap.2(as.matrix(rnaseqdata), scale = "none", col = bluered(200),key=FALSE,keysize=0.75, key.par = list(cex=0.3),
          trace = "none", density.info = "none",dendrogram = "row",Colv = FALSE, cexRow=0.9)
dev.off()
