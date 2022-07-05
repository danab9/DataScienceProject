#RNA SEQ

library(stringr)
library(DESeq2)
library(dplyr)
library(apeglm)
library(readxl)
library(pheatmap)
#library(RColorBrewer)
library(pals)

source('~/DataScienceProject/RNA-seq/differentialexpression_function.R')

# import data
rnaseqdata <- read.table("~/Dokumente/GSE153873_summary_count.star.txt", sep = "\t",header = TRUE,row.names = 1)
colnames(rnaseqdata) <- paste0(data.frame(str_split(colnames(rnaseqdata),'[.]'))[3,],data.frame(str_split(colnames(rnaseqdata),'[.]'))[1,])
#rnaseq_ercc <- read.table("RNAseq/GSE153873_summary_count.ercc.txt", sep = "\t")
dim(rnaseqdata)
# 27135 different genes and 30 patients
#rnaseqdata <- log2(rnaseqdata+1)

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



upregulated_list <- 2^(upregulated$log2FoldChange) 
names(upregulated_list)<- as.character(upregulated$genes)
upregulated_list <- sort(upregulated_list, decreasing=TRUE)

downregulated_list <- 2^(downregulated$log2FoldChange) 
names(downregulated_list)<- as.character(downregulated$genes)
downregulated_list <- sort(downregulated_list, decreasing=TRUE)


library(enrichplot)
library(DOSE)
library(org.Hs.eg.db)
library(clusterProfiler)
library(ggplot2)
#data(geneList)
up_BPs <- gseGO( upregulated_list, 
                 keyType = "SYMBOL",
                 OrgDb = org.Hs.eg.db, 
                 ont = "BP",
                 pvalueCutoff = 0.05,
                 nPermSimple = 10000 )

up_CCs <- gseGO( upregulated_list,
                 keyType = "SYMBOL",
                 OrgDb = org.Hs.eg.db, 
                 ont = "CC",
                 pvalueCutoff = 0.05,
                 nPermSimple = 10000 )

up_MFs <- gseGO( upregulated_list, 
                 keyType = "SYMBOL",
                 OrgDb = org.Hs.eg.db, 
                 ont = "MF",
                 pvalueCutoff = 0.05,
                 nPermSimple = 10000 )




down_BPs <- gseGO( downregulated_list, 
                   keyType = "SYMBOL",
                 OrgDb = org.Hs.eg.db, 
                 ont = "BP",
                 pvalueCutoff = 0.05,
                 nPermSimple = 10000 )

down_CCs <- gseGO( downregulated_list, 
                   keyType = "SYMBOL",
                 OrgDb = org.Hs.eg.db, 
                 ont = "CC",
                 pvalueCutoff = 0.05,
                 nPermSimple = 10000 )

down_MFs <- gseGO( downregulated_list, 
                   keyType = "SYMBOL",
                 OrgDb = org.Hs.eg.db, 
                 ont = "MF",
                 pvalueCutoff = 0.05,
                 nPermSimple = 10000 )





# Visualization 

# enriched GOS 
dotplot(up_CCs, showCategory=20) + ggtitle("Cellular Component upregulated GOs")

#similiar to dotplot, but GOs are clustered 
treeplot(pairwise_termsim(up_BPs), showCategory = 20, hclust_method = "average")


# distribution of log2 fold changes of enriched GOs
ridgeplot(down_MFs, showCategory=15) + xlim(c(1,2))


## network plots with genes , not working good, too many genes
## convert gene ID to Symbol
edox <- setReadable(up_MFs, 'org.Hs.eg.db', 'ENTREZID')
p1 <- cnetplot(edox, foldChange= upregulated_list)
## categorySize can be scaled by 'pvalue' or 'geneNum'
p2 <- cnetplot(edox, categorySize="pvalue", foldChange=upregulated_list)
p3 <- cnetplot(edox, foldChange=upregulated_list, circular = TRUE, colorEdge = TRUE, max.overlaps= 500) 
cowplot::plot_grid(p1, p2, p3, ncol=3, labels=LETTERS[1:3], rel_widths=c(.8, .8, 1.2))



## KEGG enrichment 

ids<-bitr(downregulated$genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb= org.Hs.eg.db)
dedup_ids = ids[!duplicated(ids[c("SYMBOL")]),]

df2 = right_join(downregulated, dedup_ids, by= c("genes" = "SYMBOL"))

kegg_gene_list <- (df2$log2FoldChange)

# Name vector with ENTREZ ids
names(kegg_gene_list) <- as.character(df2$ENTREZID)

# omit any NA values 
kegg_gene_list<-na.omit(kegg_gene_list)

# sort the list in decreasing order (required for clusterProfiler)
kegg_gene_list = sort(kegg_gene_list, decreasing = TRUE)

up_kk <- gseKEGG(geneList = kegg_gene_list,
               organism     = "hsa",
               nPermSimple        = 10000,
               minGSSize    = 1,
               maxGSSize    = 800,
               pvalueCutoff = 0.05,
               #pAdjustMethod = "BH",
               keyType       = "ncbi-geneid")

down.kk <- enrichPathway(kegg_gene_list, )


               
dotplot(up_kk, showCategory = 10, title = "Enriched Pathways" , split=".sign") + facet_grid(.~.sign)

