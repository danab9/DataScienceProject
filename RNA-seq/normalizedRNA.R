library(stringr)
library(DESeq2)
library(EDASeq)
library(biomaRt)
library(org.Hs.eg.db)
library(corrplot)

# read RNA data
rnaseqdata <- read.table("RNAseq/GSE153873_summary_count.star.txt", sep = "\t",header = TRUE,row.names = 1)
colnames(rnaseqdata) <- paste0(data.frame(str_split(colnames(rnaseqdata),'[.]'))[3,],data.frame(str_split(colnames(rnaseqdata),'[.]'))[1,])

# get ensembl id -> ID conversion from hgnc symbols to ensembl ID
mart <- useMart(dataset="hsapiens_gene_ensembl","ensembl")
ensemblnew <- getBM(filters= "hgnc_symbol", 
                    attributes=c( 'ensembl_gene_id','hgnc_symbol'),
                    values= rownames(rnaseqdata),
                    mart= mart)
# not every hgnc symbol had an ensembl ID
ensemblnew <- ensemblnew[!is.na(ensemblnew$ensembl_gene_id),]

# get gene length
genelength <- getGeneLengthAndGCContent(id=ensemblnew$ensembl_gene_id,'hsa')

# DESeqData
condition <- factor(colnames(rnaseqdata))
deseqrna <- DESeqDataSetFromMatrix(rnaseqdata[which(!is.na(ensemblnew$ensembl_gene_id)),], DataFrame(condition), ~ condition)
mcols(deseqrna)$basepairs <- genelength[,1]

# Counts normalized per kilobase
normalizedRNA <- log2(fpkm(deseqrna, robust = TRUE)+1)

# save the normalized RNA
write.table(normalizedRNA, 'normalizedRNA.txt', quote=FALSE, append = FALSE, sep = " ", dec = ".",
            row.names = TRUE, col.names = TRUE)

##EXTRA plots for visualization
## bp of non-normalized

boxplot(log2(assay(deseqrna)+1), notch=TRUE,
        main = "Non-normalized read counts",
        ylab="log2(read counts+1)", cex = .6)
# description: For each patient the non-normalized read counts. The median is around 6.
## bp of normalized values
boxplot(normalizedRNA, notch=TRUE,
        main = "Normalized read counts per kb",
        ylab="normalized read counts", cex = .6)

# description: For each patient the normalized read counts. The median is around 2. 

library(ggplot2)
# boxplot function with wilcox test
boxplotfunction <- function(ndata,n=2){
  AD <- rowMeans(dplyr::select(ndata, starts_with("AD")))
  old <- rowMeans(dplyr::select(ndata, starts_with("Old")))
  if(n==3){
    young <- rowMeans(dplyr::select(ndata, starts_with("Young")))
    condition <- c(rep('young',length(young)),rep('old',length(old)),rep('AD',length(AD)))
    data <- data.frame(class=rep(c('young','old','AD'),c(length(young),length(old),length(AD))),normalized_RNA=c(young,old,AD))
    test <- kruskal.test(data$normalized_RNA~data$class)
    name <- "Kruskal-Wallis: "
  }else{
    condition <- c(rep('old',length(old)),rep('AD',length(AD)))
    data <- data.frame(class=rep(c('old','AD'),c(length(old),length(AD))),normalized_RNA=c(old,AD))
    test <- wilcox.test(data$normalized_RNA ~ data$class)
    name <- "Wilcox: "
  }
  annotations <- data.frame(
    xpos = -Inf,
    ypos = Inf,
    annotateText = paste(name,signif(test$p.value,digits=3)),
    hjustvar = 0 ,
    vjustvar = 1) #<- adjust
  p <- ggplot(data, aes(x=class, y=normalized_RNA),size=10) +
    geom_boxplot()+
    geom_label(data=annotations,aes(x=xpos,y=ypos,hjust=hjustvar,vjust=vjustvar,label=annotateText),size=6)+
    theme_grey(base_size = 22)
  p
  
}

normalizedRNA <- read.table("normalizedRNA.txt",header=TRUE,sep = ",",row.names = 1)
#upregulatedAll <- read.table("upregulatedAll.txt",header=TRUE,row.names = 1)
#downregulatedAll <- read.table("downregulatedAll.txt",header=TRUE,row.names = 1)
upregulated <- read.table("upregulated.txt",header=TRUE,row.names = 1)
downregulated <- read.table("downregulated.txt",header=TRUE,row.names = 1)


# new <- normalizedRNA[upregulatedAll$genes,]
# new <- new[!is.na(new),]
# boxplotfunction(new,n=3)
# 
# new <- normalizedRNA[downregulatedAll$genes,]
# new <- new[!is.na(new),]
# boxplotfunction(new,n=3)

#Boxplot of normalized read counts from upregulated/downregulated genes: old vs AD
new <- normalizedRNA[upregulated$genes,]
new <- new[!is.na(new),]
boxplotfunction(new,n=2)

new <- normalizedRNA[downregulated$genes,]
new <- new[!is.na(new),]
boxplotfunction(new,n=2)

# description: Boxplot of up/downregulated genes in old and AD group
# in addition a wilcox-test was calculated
