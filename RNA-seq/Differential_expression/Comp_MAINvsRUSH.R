library(stringr)
source('differentialexpression_function.R')

#RUSH data
diffRUSH <- read.table('diff_genes_maleRush.txt',header=TRUE,sep=',')
diffRUSHall <- read.table("diff_genes_allRUSH.txt",header = TRUE,sep=",")

# Main data up and downregulated genes
upregulated <- read.table("upregulated.txt",header=TRUE)
downregulated <- read.table("downregulated.txt",header=TRUE)

# male RUSH
which(downregulated$genes%in%diffRUSH$hgnc_symbol[diffRUSH$hgnc_symbol!=""])
which(upregulated$genes%in%diffRUSH$hgnc_symbol[diffRUSH$hgnc_symbol!=""])
# the differential expressed genes from the RUSH data are not significant up or downregulated
# in the main data

# all RUSH
which(downregulated$genes%in%diffRUSHall$hgnc_symbol[diffRUSH$hgnc_symbol!=""])
which(upregulated$genes%in%diffRUSHall$hgnc_symbol[diffRUSH$hgnc_symbol!=""])
# same as before

# redoing the differential expression to check if the log2-fold character is the same
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
plotdf <- plotdf[order(plotdf$genes),]

d1 <- plotdf[which(plotdf$genes%in%diffRUSH$hgnc_symbol[diffRUSH$hgnc_symbol!=""]),]
d1_1 <- diffRUSH[diffRUSH$hgnc_symbol%in%d1$genes,c(1,2,5)]
cbind(d1,d1_1[order(d1_1$hgnc_symbol),])

d2 <- plotdf[which(plotdf$genes%in%diffRUSHall$hgnc_symbol[diffRUSHall$hgnc_symbol!=""]),]
d2_2 <- diffRUSHall[diffRUSHall$hgnc_symbol%in%d2$genes,c(1,2,5)]
cbind(d2,d2_2[order(d2_2$hgnc_symbol),])
# sometimes the log2-fold is in some cases both negativ/positiv