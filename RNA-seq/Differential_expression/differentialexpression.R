#RNA SEQ
library(stringr)
library(DESeq2)
library(dplyr)
library(apeglm)
library(readxl)
library(pheatmap)
library(pals)
library(PerformanceAnalytics)
source('differentialexpression_function.R')


# import data
rnaseqdata <- read.table("RNAseq/GSE153873_summary_count.star.txt", sep = "\t",header = TRUE,row.names = 1)
colnames(rnaseqdata) <- paste0(data.frame(str_split(colnames(rnaseqdata),'[.]'))[3,],data.frame(str_split(colnames(rnaseqdata),'[.]'))[1,])
#rnaseq_ercc <- read.table("RNAseq/GSE153873_summary_count.ercc.txt", sep = "\t")
dim(rnaseqdata)
# 27135 different genes and 30 patients
any(is.na(rnaseqdata))

# Groups
AD <- dplyr::select(rnaseqdata, starts_with("AD"))
young <- dplyr::select(rnaseqdata, starts_with("Young"))
old <- dplyr::select(rnaseqdata, starts_with("Old"))

# new order of patient
rnaseqdata <- data.frame(AD,old,young)


# First task: Comparing gene expression between AD and old samples

# DIFFERENTIAL EXPRESSION
res <- diffexpression(old,AD,alpha=0.05,tidy=FALSE)
summary(res)
dds <- diffexpression(old,AD,alpha=0.05,tidy=FALSE,result=FALSE)
#condition <- factor(c(rep("old",ncol(old)),rep("AD",ncol(AD))))
#421 genes with significant upregulation in AD, while 434
# had significant downregulation -> 855 significant up/downregulated

resOrdered <- res[order(res$pvalue),]
#Yekutieli control 

pvalues <- p.adjust(res$pvalue, method = "BY")
sum(pvalues<0.05,na.rm = T)


plotdf<- data_frame(log2FoldChange=res$log2FoldChange, padj= res$padj, Yekutielipval=pvalues,pval=res$pvalue, genes=rownames(res))
plotdf<- plotdf %>%
  mutate(gene_type = case_when(log2FoldChange > 0 & padj < 0.05 ~ "upregulated",
                               log2FoldChange < 0 & padj < 0.05 ~ "downregulated",
                               TRUE ~ "not significant"))   
print(table(plotdf$gene_type))

#up and downregulated genes, ordered by p_adjust
upregulated <- plotdf[plotdf$gene_type=="upregulated",]
upregulated <- upregulated[order(upregulated$padj),]
downregulated <- plotdf[plotdf$gene_type=="downregulated",]
downregulated <- downregulated[order(downregulated$padj),]



# Comparison to the paper's results
#434 downregulated, 421 upregulated
sheets <- excel_sheets("RNAseq/GSE153873_AD.vs.Old_diff.genes.xlsx")
resultupreg <- read_excel("RNAseq/GSE153873_AD.vs.Old_diff.genes.xlsx",sheet = sheets[1])
resultdownreg <- read_excel("RNAseq/GSE153873_AD.vs.Old_diff.genes.xlsx",sheet = sheets[2])


table(plotdf[plotdf$genes %in% unlist(resultupreg[,1]),"gene_type"])
table(plotdf[plotdf$genes %in% unlist(resultdownreg[,1]),"gene_type"])


# Save up-/downregulated count matrix and up-/downregulated genes 
write.table(as.factor(unlist(upregulated$genes)), 'upregulatedgenes.txt',quote = FALSE, append = FALSE, sep = " ", dec = ".",
            row.names = FALSE, col.names = FALSE)

write.table(rnaseqdata[rownames(rnaseqdata)%in%upregulated$genes,], 'upregulated_count.txt', quote=FALSE, append = FALSE, sep = " ", dec = ".",
            row.names = TRUE, col.names = TRUE)


write.table(upregulated, 'upregulated_results.txt', quote=FALSE, append = FALSE, sep = " ", dec = ".",
            row.names = TRUE, col.names = TRUE)

write.table(downregulated$genes, 'downregulationgenes.txt',quote = FALSE, append = FALSE, sep = " ", dec = ".",
            row.names = FALSE, col.names = FALSE)

write.table(downregulated, 'downregulated_results.txt', quote = FALSE, append = FALSE, sep = " ", dec = ".",
            row.names = TRUE, col.names = TRUE)

#three main genes from the paper
imp <- c( 'EP300', 'CREBBP','TRRAP')



################################################################################
## VISUALIZAION


# Plot counts for a single gene -> smallest pvalue
plotCounts(dds, gene=which.min(res$padj), intgroup="condition")
# MA plot
png(file = "MAplot.png", height = 400,width=400)
plotMA(res, ylim=c(-3,3),colSig='red',colLine='red',alpha=0.05,cex=0.4,xlab='Mean expression',ylab='log2 (AD/old)')
dev.off()
# description: MA plot of differential results (AD vs old)

#Alternative shrinkage estimators
resNorm <- lfcShrink(dds, coef=2, type="normal")
plotMA(resNorm, ylim=c(-1.5,1.5),main='Normal')



#Heatmap of the count matrix + cluster
upregulated <- upregulated[order(upregulated$log2FoldChange,decreasing = TRUE),]
ntd <- assay(normTransform(dds))
ntd <- assay(dds)
select <- order(rowMeans(counts(dds,normalized=TRUE)),
                decreasing=TRUE)[1:60]

pheatmap(ntd[which(rownames(ntd) %in% upregulated$genes[1:60]),],cluster_cols = FALSE)

library("gplots")
x11()
heatmap.2(log2(1+ntd[which(rownames(ntd) %in% upregulated$genes[1:60]),]), scale = "none", col = bluered(200), 
          trace = "none", density.info = "none",dendrogram = "row")
##Principal component plot of the samples
#Related to the distance matrix is the PCA plot, which shows the samples in the
#2D plane spanned by their first two principal components. This type of plot is
#useful for visualizing the overall effect of experimental covariates and batch effects.
vsd <- vst(dds, blind=FALSE)
vsd$condition <- factor(c(rep("old",ncol(old)),rep("AD",ncol(AD))))
levels(vsd$condition) <- c('old','AD')
png(file = "PCARNAseq.png", height = 400,width=600)
plotPCA(vsd, intgroup=c("condition"))
dev.off()
# description: PCA plot AD vs old 


# other differential expression
res <- diffexpression(young,AD,alpha=0.05,tidy=FALSE)
summary(res)

res <- diffexpression(young,old,alpha=0.05,tidy=FALSE)
summary(res)

res <- diffexpression(data.frame(young,old),AD,alpha=0.05,tidy=FALSE)
summary(res)

plotdf<- data_frame(log2FoldChange=res$log2FoldChange, padj= res$padj,pval=res$pvalue, genes=rownames(res))
plotdf<- plotdf %>%
  mutate(gene_type = case_when(log2FoldChange > 0 & padj < 0.05 ~ "upregulated",
                               log2FoldChange < 0 & padj < 0.05 ~ "downregulated",
                               TRUE ~ "not significant"))   
print(table(plotdf$gene_type))

#up and downregulated genes, ordered by p_adjust
upregulatedAll <- plotdf[plotdf$gene_type=="upregulated",]
upregulatedAll <- upregulatedAll[order(upregulatedAll$padj),]
downregulatedAll <- plotdf[plotdf$gene_type=="downregulated",]
downregulatedAll <- downregulatedAll[order(downregulatedAll$padj),]

write.table(upregulatedAll, 'upregulatedAll.txt', quote=FALSE, append = FALSE, sep = " ", dec = ".",
            row.names = TRUE, col.names = TRUE)

write.table(downregulatedAll, 'downregulatedAll.txt',quote = FALSE, append = FALSE, sep = " ", dec = ".",
            row.names = TRUE, col.names = TRUE)


