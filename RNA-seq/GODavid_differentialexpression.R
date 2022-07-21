library(stringr)
library(DESeq2)
library(dplyr)
library(apeglm)
library(readxl)
library(pheatmap)
library(mixOmics)
library(biomaRt)
library(gplots)

source('differentialexpression_function.R')

# import RNA-seq data
rnaseqdata <- read.table("RNAseq/GSE153873_summary_count.star.txt", sep = "\t",header = TRUE,row.names = 1)
colnames(rnaseqdata) <- paste0(data.frame(str_split(colnames(rnaseqdata),'[.]'))[3,],data.frame(str_split(colnames(rnaseqdata),'[.]'))[1,])

# import genes from GO enrichment "Regulation of transcription"
genes_regtrans <- read.table("RNAseq/GOenrichment/reg_transcription_own.txt", sep = "\t",header = TRUE)
colnames(genes_regtrans) <- c("Gene Name","Species")

# save the Gene names for STRING analysis
write.table(rownames(genes_regtrans), 'go_genes.txt', quote=FALSE, append = FALSE, sep = " ", dec = ".",
              row.names = FALSE, col.names = FALSE)

# load the genes from GO term "regulation of transcription" from the paper
genes_regtrans_paper <- read.table("RNAseq/GOenrichment/reg_transcription_paper.txt", sep = "\t",header = TRUE)
colnames(genes_regtrans_paper) <- c("Gene Name","Species")
# save only the gene names
write.table(rownames(genes_regtrans_paper), 'go_genesPaper.txt', quote=FALSE, append = FALSE, sep = " ", dec = ".",
            row.names = FALSE, col.names = FALSE)

# overlap of genes in "regulation of transcription" own and paper results
length(intersect(rownames(genes_regtrans),rownames(genes_regtrans_paper)))

# differential expression with the "regulation of transcription" genes

ownres <- unlist(read.table( 'go_genes.txt'))
#ownres <- ownres[!ownres%in%c('LINC-PINT','SGF29','TFDP2','ZNF137P')]

rnaseqdata <- rnaseqdata[rownames(rnaseqdata)%in%ownres,]
AD <-  dplyr::select(rnaseqdata, starts_with("AD"))
young <-  dplyr::select(rnaseqdata, starts_with("Young"))
old <-  dplyr::select(rnaseqdata, starts_with("Old"))

rnaseqdata <- data.frame(young,old,AD)


# Comparing gene expression between young AD and old samples

dds <- diffexpression(young,old,AD,result=FALSE)
dds


#Heatmap of the count matrix + cluster
# normalized counts transformation
ntd <- assay(normTransform(dds))

png(file = "RT_heatmap.png", height = 1800) 
heatmap.2(ntd, scale = "none", col = bluered(200),keysize=0.75, key.par = list(cex=0.3),key=FALSE,
          trace = "none", density.info = "none",dendrogram = "row",Colv = FALSE, cexRow=0.9)
dev.off()
# Description: Heatmap of normalized counts with DESeq from differential expressed genes. 
# These genes are within the GO term "regulation of transcription". Only the rows with the genes are clustered (same as in the paper).
# X-axis: young, old and AD group, Y-axis: different genes within "regulation of transcription"


# Heatmap of normalized RNA with GO term genes
RNAnorm <- read.table("RNAseq/normalizedRNA.txt")
RNAnorm <- RNAnorm[rownames(RNAnorm)%in%ownres,]
AD <-  dplyr::select(RNAnorm, starts_with("AD"))
young <-  dplyr::select(RNAnorm, starts_with("Young"))
old <-  dplyr::select(RNAnorm, starts_with("Old"))

RNAnorm <- data.frame(young,old,AD)

png(file = "RT_heatmapNew.png", height = 1800)
heatmap.2(as.matrix(RNAnorm), scale = "none", col = bluered(200),key=FALSE,keysize=0.75, key.par = list(cex=0.3),
          trace = "none", density.info = "none",dendrogram = "row",Colv = FALSE, cexRow=0.9)
dev.off()
# Description: Heatmap of own normalized RNA (see normalizedRNA.R). Used genes within the GO term "regulation of transcription". 
# Only the rows with the genes are clustered (same as in the paper).
# X-axis: young, old and AD group, Y-axis: different genes within "regulation of transcription"

#heatmap.2(as.matrix(RNAnorm), scale = "none", col = bluered(200),key=FALSE,keysize=0.75, key.par = list(cex=0.3),
#          trace = "none", density.info = "none", cexRow=0.9)
