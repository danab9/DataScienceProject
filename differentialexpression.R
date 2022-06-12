#RNA SEQ


library(stringr)
library(DESeq2)
library(dplyr)
library(apeglm)
library(readxl)
library(pheatmap)


# import data
rnaseqdata <- read.table("RNAseq/GSE153873_summary_count.star.txt", sep = "\t",header = TRUE,row.names = 1)
colnames(rnaseqdata) <- paste0(data.frame(str_split(colnames(rnaseqdata),'[.]'))[3,],data.frame(str_split(colnames(rnaseqdata),'[.]'))[1,])
#rnaseq_ercc <- read.table("RNAseq/GSE153873_summary_count.ercc.txt", sep = "\t")
dim(rnaseqdata)
# 27135 different genes and 30 patients


# Groups
AD <- select(rnaseqdata, starts_with("AD"))
young <- select(rnaseqdata, starts_with("Young"))
old <- select(rnaseqdata, starts_with("Old"))

# new order of patient
rnaseqdata <- data.frame(AD,old,young)


# First task: Comparing gene expression between AD and old samples
# http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html
condition <- factor(c(rep("Old",ncol(old)),rep("AD",ncol(AD))))
dds <- DESeqDataSetFromMatrix(countData = data.frame(old,AD),
                              DataFrame(condition), ~ condition)
dds <- DESeq(dds)
res <- results(dds,alpha=0.05,tidy=TRUE) # to obtain also the genes


#421 genes with significant upregulation in AD, while 434
# had significant downregulation -> 855 significant up/downregulated
resOrdered <- res[order(res$pvalue),]
summary(res)

pvalues <- p.adjust(res$pvalue, method = "BH")
sum(pvalues<0.05,na.rm = T)


# MA plot
plotMA(res, ylim=c(-4,4))

resNorm <- lfcShrink(dds, coef=2, type="normal")
plotMA(resNorm, ylim=c(-1.5,1.5),main='Normal')


#434 downregulated, 421 upregulated
plotdf<- data_frame(log2FoldChange=res$log2FoldChange, padj= res$padj,genes=res$row)

plotdf<- plotdf %>%
  mutate(gene_type = case_when(log2FoldChange >= 0.6999 & padj < 0.05 ~ "upregulated",
                               log2FoldChange <= -0.668 & padj < 0.05 ~ "downregulated",
                               TRUE ~ "not significant"))   
# OR
plotdf<- data_frame(log2FoldChange=res$log2FoldChange, padj= res$padj)
plotdf<- plotdf %>%
  mutate(gene_type = case_when(log2FoldChange > 0.661 & pvalues < 0.05 ~ "upregulated",
                               log2FoldChange < -0.615 & pvalues < 0.05 ~ "downregulated",
                               TRUE ~ "not significant"))   

print(table(plotdf$gene_type))

upregulated <- plotdf[plotdf$gene_type=="upregulated",]
downregulated <- plotdf[plotdf$gene_type=="downregulated",]

# Comparison to the paper's results
sheets <- excel_sheets("RNAseq/GSE153873_AD.vs.Old_diff.genes.xlsx")
x <- lapply(sheets, function(X) readxl::read_excel(filename, sheet = X))
resultupreg <- read_excel("RNAseq/GSE153873_AD.vs.Old_diff.genes.xlsx",sheet = sheets[1])
resultdownreg <- read_excel("RNAseq/GSE153873_AD.vs.Old_diff.genes.xlsx",sheet = sheets[2])




# Plot counts for a single gene -> smallest pvalue
plotCounts(dds, gene=which.min(res$padj), intgroup="condition")

#Heatmap of the count matrix + cluster
ntd <- normTransform(dds)
select <- order(rowMeans(counts(dds,normalized=TRUE)),
                decreasing=TRUE)[1:15]
pheatmap(assay(ntd)[select,])

#Principal component plot of the samples

plotPCA(vsd, intgroup=c("condition"))