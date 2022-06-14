library(DESeq2)
library(readr)
library(stringr)
library(dplyr)
rush_data<- read.csv("/home/rosa/rush_data/all_counts.csv",  header=TRUE)
rush_data<- rush_data[,2:ncol(rush_data)]
rush_data<- dplyr::select(rush_data, -c(ensembl_gene_id, gene_biotype))
rownames(rush_data)<- rush_data$hgnc_symbol


rnaseqdata <- read.table("/home/rosa/Dokumente/GSE153873_summary_count.star.txt", sep = "\t",header = TRUE,row.names = 1)
colnames(rnaseqdata) <- paste0(data.frame(str_split(colnames(rnaseqdata),'[.]'))[3,],data.frame(str_split(colnames(rnaseqdata),'[.]'))[1,])
#rnaseq_ercc <- read.table("RNAseq/GSE153873_summary_count.ercc.txt", sep = "\t")
dim(rnaseqdata)
# 27135 different genes and 30 patients


# Groups
AD <- dplyr::select(rnaseqdata, starts_with("AD"))
young <- dplyr::select(rnaseqdata, starts_with("Young"))
old <- dplyr::select(rnaseqdata, starts_with("Old"))
female<- dplyr::select(rush_data, -hgnc_symbol)

old<- dplyr::select(old, -hgnc_symbol)


# differential gene expression: control vs female
# 23589 common measured genes 

control_female<- inner_join(old, rush_data,  by="hgnc_symbol")
rownames(control_female)<- control_female$hgnc_symbol
control_female<- dplyr::select(control_female, -hgnc_symbol)

AD$hgnc_symbol<- rownames(AD)
male_female<- inner_join(AD, rush_data,  by="hgnc_symbol")
rownames(male_female)<- male_female$hgnc_symbol
male_female<- dplyr::select(male_female, -hgnc_symbol)
AD<-  dplyr::select(AD, -hgnc_symbol)


execute_dge<- function(matrix, samplesize1, samplesize2){
  condition <- factor(c(rep("G1",samplesize1),rep("G2",samplesize2)))
  
  dds_cf <- DESeqDataSetFromMatrix(countData = matrix,
                                   DataFrame(condition), ~ condition)
  dds_cf <- DESeq(dds_cf)
  return (dds_cf)
  
}

dds<- execute_dge(male_female, ncol(AD), ncol(female))
res<- results(dds,alpha= 0.05)

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
