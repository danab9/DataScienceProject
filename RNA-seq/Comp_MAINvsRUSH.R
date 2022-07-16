library(stringr)
source('differentialexpression_function.R')

diffRUSH <- read.table('diff_genes_maleRush.txt',header=TRUE,sep=',')
diffRUSHall <- read.table("diff_genes_allRUSH.txt",header = TRUE,sep=",")

upregulated <- read.table("upregulated.txt",header=TRUE)
downregulated <- read.table("downregulated.txt",header=TRUE)

which(downregulated$genes%in%diffRUSH$hgnc_symbol[diffRUSH$hgnc_symbol!=""])

rnaseqdata <- read.table("RNAseq/GSE153873_summary_count.star.txt", sep = "\t",header = TRUE,row.names = 1)
colnames(rnaseqdata) <- paste0(data.frame(str_split(colnames(rnaseqdata),'[.]'))[3,],data.frame(str_split(colnames(rnaseqdata),'[.]'))[1,])
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
plotdf<- data.frame(log2FoldChange=res$log2FoldChange, padj= res$padj, genes=rownames(res))
plotdf[which(plotdf$genes%in%diffRUSH$hgnc_symbol[diffRUSH$hgnc_symbol!=""]),]
plotdf[which(plotdf$genes%in%diffRUSHall$hgnc_symbol[diffRUSHall$hgnc_symbol!=""]),]
