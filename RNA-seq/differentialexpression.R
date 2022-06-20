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
#rnaseqdata <- log2(rnaseqdata+1)

# Groups
AD <- select(rnaseqdata, starts_with("AD"))
young <- select(rnaseqdata, starts_with("Young"))
old <- select(rnaseqdata, starts_with("Old"))

# new order of patient
rnaseqdata <- data.frame(AD,old,young)


# First task: Comparing gene expression between AD and old samples
# http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html

diffexpression <- function(group1,group2,alpha,tidy=TRUE){
  condition <- factor(c(rep("G1",ncol(group1)),rep("G2",ncol(group2))))
  dds <- DESeqDataSetFromMatrix(countData = data.frame(group1,group2),
                              DataFrame(condition), ~ condition)
  dds <- DESeq(dds)
  
  return(results(dds,alpha=alpha,tidy=tidy)) # to obtain also the genes
}

res <- diffexpression(old,AD,alpha=0.05,tidy=FALSE)
summary(res)

#421 genes with significant upregulation in AD, while 434
# had significant downregulation -> 855 significant up/downregulated

resOrdered <- res[order(res$pvalue),]
#Yekutieli control 

pvalues <- p.adjust(res$pvalue, method = "BY")
sum(pvalues<0.05,na.rm = T)


# MA plot
plotMA(res, ylim=c(-4,4))

resNorm <- lfcShrink(res, coef=2, type="normal")
plotMA(resNorm, ylim=c(-1.5,1.5),main='Normal')


#434 downregulated, 421 upregulated
# plotdf<- data_frame(log2FoldChange=res$log2FoldChange, padj= res$padj,genes=rownames(res))
# 
# # OR
# plotdf<- data_frame(log2FoldChange=res$log2FoldChange, padj= res$padj)
# plotdf<- plotdf %>%
#   mutate(gene_type = case_when(log2FoldChange > 0.661 & pvalues < 0.05 ~ "upregulated",
#                                log2FoldChange < -0.615 & pvalues < 0.05 ~ "downregulated",
#                                TRUE ~ "not significant"))   
# 
# plotdf<- plotdf %>%
#   mutate(gene_type = case_when(log2FoldChange >= 0.681 & padj < 0.05 ~ "upregulated",
#                                log2FoldChange <= -0.6869 & padj < 0.05 ~ "downregulated",
#                                TRUE ~ "not significant"))   

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


# Save up-/downregulated count matrix and up-/downregulated genes 
write.table(as.factor(unlist(upregulated$genes)), 'upregulatedgenes.txt',quote = FALSE, append = FALSE, sep = " ", dec = ".",
            row.names = FALSE, col.names = FALSE)

write.table(rnaseqdata[rownames(rnaseqdata)%in%upregulated$genes,], 'upregulated_count.txt', quote=FALSE, append = FALSE, sep = " ", dec = ".",
            row.names = TRUE, col.names = TRUE)

write.table(downregulated$genes, 'downregulationgenes.txt',quote = FALSE, append = FALSE, sep = " ", dec = ".",
            row.names = FALSE, col.names = FALSE)

write.table(rnaseqdata[rownames(rnaseqdata)%in%upregulated$genes,], 'upregulated_count.txt', quote = FALSE, append = FALSE, sep = " ", dec = ".",
            row.names = TRUE, col.names = TRUE)


# Comparison to the paper's results
#434 downregulated, 421 upregulated
sheets <- excel_sheets("RNAseq/GSE153873_AD.vs.Old_diff.genes.xlsx")
resultupreg <- read_excel("RNAseq/GSE153873_AD.vs.Old_diff.genes.xlsx",sheet = sheets[1])
resultdownreg <- read_excel("RNAseq/GSE153873_AD.vs.Old_diff.genes.xlsx",sheet = sheets[2])

# overlap of up-/downregulated genes in own results and papers results  
length(intersect(upregulated$genes,unlist(resultupreg[,1])))
#314 of 421
length(intersect(downregulated$genes,unlist(resultdownreg[,1])))
#344 of 434

table(plotdf[plotdf$genes %in% unlist(resultupreg[,1]),"gene_type"])
table(plotdf[plotdf$genes %in% unlist(resultdownreg[,1]),"gene_type"])

# a lot of overlaps but also some up/downregulated genes but also a lot of genes missing
# in the 50 most significant downregulated genes: overlap of 49 genes
# in the 50 most significant upregulated genes: overlap of 45 genes

#three main genes from the paper
imp <- c( 'EP300', 'CREBBP','TRRAP')

# Plot counts for a single gene -> smallest pvalue
plotCounts(dds, gene=which.min(res$padj), intgroup="condition")

#Heatmap of the count matrix + cluster
ntd <- normTransform(dds)
select <- order(rowMeans(counts(dds,normalized=TRUE)),
                decreasing=TRUE)[1:60]
pheatmap(assay(ntd)[select,])

#Principal component plot of the samples

plotPCA(vsd, intgroup=c("condition"))
