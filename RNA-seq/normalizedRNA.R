library(stringr)
library(DESeq2)
library(EDASeq)
library('org.Hs.eg.db')

# data
rnaseqdata <- read.table("RNAseq/GSE153873_summary_count.star.txt", sep = "\t",header = TRUE,row.names = 1)
colnames(rnaseqdata) <- paste0(data.frame(str_split(colnames(rnaseqdata),'[.]'))[3,],data.frame(str_split(colnames(rnaseqdata),'[.]'))[1,])

# get entrez id
ensembl <- mapIds(org.Hs.eg.db, rownames(rnaseqdata), 'ENSEMBL', 'SYMBOL')
ensemblnew <- ensembl[!is.na(ensembl)]
genelength <- getGeneLengthAndGCContent(id=ensemblnew,'hsa')

# DESeqData
condition <- factor(colnames(rnaseqdata))
deseqrna <- DESeqDataSetFromMatrix(rnaseqdata[which(!is.na(ensembl)),], DataFrame(condition), ~ condition)
mcols(deseqrna)$basepairs <- genelength[,1]

# Counts normalized per kilobase
normalizedRNAseq <- log2(fpkm(deseqrna, robust = TRUE)+1)

#save normalized RNA in a table
write.table(normalizedRNAseq, 'normalizedRNA.txt', quote=FALSE, append = FALSE, sep = " ", dec = ".",
            row.names = TRUE, col.names = TRUE)

##EXTRA plots for visualization
## bp of non-normalized
boxplot(log2(deseqrna+1), notch=TRUE,
        main = "Non-normalized read counts",
        ylab="log2(read counts)", cex = .6)
## bp of normalized values
boxplot(log2(normalizedRNAseq +1), notch=TRUE,
        main = "Normalized read counts per kb",
        ylab="log2(read counts)", cex = .6)


