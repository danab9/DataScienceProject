## Differential Gene Expression Analysis For Rush 
library(tximport)
library(dplyr)
library(ggplot2)
## importing files 
dir.ADM<- "/home/rosa/rush_data/ADmale"
dir.ADF<- "/home/rosa/rush_data/ADfemale"
dir.CF <- "/home/rosa/rush_data/ControlFemale"
dir.CM <- "/home/rosa/rush_data/ControlMale"

overview.ADM <- read.csv("/home/rosa/DataScienceProject/RushData/ADmale_overview.csv", header=TRUE)
overview.ADF <- read.csv("/home/rosa/DataScienceProject/RushData/ADfemale_overview.csv", header=TRUE)
overview.CF <- read.csv("/home/rosa/DataScienceProject/RushData/ControlFemale_overview.csv", header=TRUE)
overview.CM <- read.csv("/home/rosa/DataScienceProject/RushData/ControlMale_overview.csv", header=TRUE)


files.ADM <- file.path(dir.ADM, paste0(overview.ADM$file_name, ".tsv"))
names(files.ADM) <-overview.ADM$Accession
files.ADF <- file.path(dir.ADF, paste0(overview.ADF$file_name, ".tsv"))
names(files.ADF) <-overview.ADF$Accession
files.CF <- file.path(dir.CF, paste0(overview.CF$file_name, ".tsv"))
names(files.CF) <-overview.CF$Accession
files.CM <- file.path(dir.CM, paste0(overview.CM$file_name, ".tsv"))
names(files.CM) <-overview.CM$Accession







## age to 0: under 90, 1: over 90 
change.age.column <- function(overview, binary=FALSE){
  overview$age <- as.numeric(gsub("[A-z]","", overview$Biosample.age))
  if (binary){
    
    overview <- mutate(overview, age= age < 90 )
 
  }
  return(overview)
  }

overview.CF <- change.age.column(overview.CF,T)
overview.CM <- change.age.column(overview.CM,T)
overview.ADF <- change.age.column(overview.ADF,T)
overview.ADM <- change.age.column(overview.ADM,T)




# number of samples
number.controlsamples <- length(files.CF)+ length(files.CM)
number.ADsamples<-  length(files.ADF)+ length(files.ADM)


# combining all rna seq data to one rsem object
files.all <- c(files.CF, files.CM, files.ADF, files.ADM)
txi.rsem <- tximport(files.all, type = "rsem", txIn = FALSE, txOut = FALSE)

# pre-filtering genes with length 0 and count 0 
zero_length_and_unexpressed = (apply(txi.rsem$abundance, 1, max) == 0) &
  (apply(txi.rsem$length, 1, min) == 0)

txi.rsem$length = txi.rsem$length[!zero_length_and_unexpressed,]
txi.rsem$abundance = txi.rsem$abundance[!zero_length_and_unexpressed,]
txi.rsem$counts = txi.rsem$counts[!zero_length_and_unexpressed,]

# create col data : condition sex, age of samples

condition<- factor(c(rep("C", length(files.CF) + length(files.CM) ),
                    rep("AD",length(files.ADF)+  length(files.ADM))
  
))

sex<- factor(c(rep("F", length(files.CF)),
               rep("M", length(files.CM)), 
               rep("F",length(files.ADF)),
               rep("M", length(files.ADM))
))

age<- factor(c( overview.CF$age,  overview.CM$age, overview.ADF$age, overview.ADM$age ))


sampletable <- data.frame(condition= condition, sex=sex, age=age)


# create dds object (summarizedexperiment obj)

dds <- DESeqDataSetFromTximport(txi = txi.rsem ,colData = sampletable ,design = ~ age + condition)
#pre-filtering 
# ERCC counts, no genes
dds<- dds[1:(nrow(dds)-90),]
# low gene expression measurements across all samples
dds <- dds[ rowSums(counts(dds)) > 1, ]

dds.female <- dds[, which(sex=="F")]
dds.male<- dds[, which(sex=="M")] #  & 
dds.male <- DESeq(dds[, which(sex=="M")])
# data exploration

# variance stabilizing transformation. controlling great variance across the mean ( negative binomial data)
vsd<- vst(dds.male, blind=FALSE)
# or  regularized-logarithm transformation, takes more time
#rld <- rlog(dds, blind = FALSE)


dds <- estimateSizeFactors(dds)

create_log2_comparsion <- function(rows){
  
  df <- bind_rows(
    as_data_frame(log2(counts(dds, normalized=TRUE)[, rows]+1)) %>%
      mutate(transformation = "log2(x + 1)"),
    as_data_frame(assay(vsd)[, 1:2]) %>% mutate(transformation = "vst"),
  #  as_data_frame(assay(rld)[, 1:2]) %>% mutate(transformation = "rlog")
  )
  
  colnames(df)[1:2] <- c("x", "y")  
  
  plot<- ggplot(df, aes(x = x, y = y))+ geom_hex(bins = 80) +
  coord_fixed() + facet_grid( . ~ transformation)  
  
  return(plot)
}

plot_cf<- create_log2_comparsion(c(1,2))
plot_cm<- create_log2_comparsion(c(27,28))
plot_control<- create_log2_comparsion(c(1, 27))
plot_af <- create_log2_comparsion(c(42,43))
plot_am <- create_log2_comparsion(c(71,72))
plot_AD<- create_log2_comparsion(c(42,71))


# sample sample distance 

#PCA

pcaData <- plotPCA(vsd, intgroup= c("condition", "age"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=age, shape=condition)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) 


# sex adjustment, batch/sex removal 

assay(vsd) <- limma::removeBatchEffect(assay(vsd), vsd$sex)

pcaData <- plotPCA(vsd, intgroup= c("condition", "sex"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=condition, shape=sex)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) 


# using combatseq to correct for sex
combatseq_counts<- ComBat_seq(as.matrix(assay(dds.all)), batch = sex, group=condition )
# PCA 
pca_combatseq<- as.data.frame(prcomp(combatseq_counts)[2]$rotation)
pca_combatseq$condition<- condition 
pca_combatseq$batch<- batch
ggplot(pca_combatseq, aes(PC1, PC2, color=condition, shape = sex)) +
  geom_point(size=3) 

# still not showing distance between AD and control samples


# heatmap 

library("pheatmap")
library("RColorBrewer")
sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd$condition, vsd$sex, sep="-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pdf("heatmap_sex_removal.pdf", height = 4, width = 5)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)
dev.off()



### differential gene expression analysis

res<- results(dds.male,alpha=0.05)

resOrdered <- res[order(res$pvalue),]
summary(res)

pvalues <- p.adjust(res$pvalue, method = "BH")
sum(pvalues<0.05,na.rm = T)

# only 9  diff. expr. gene detected !!!

# MA plot
plotMA(res, ylim=c(-4,4))

resNorm <- lfcShrink(dds.male, coef=2, type="normal")
plotMA(resNorm, ylim=c(-1.5,1.5),main='Normal')



plotdf<- data_frame(log2FoldChange=res$log2FoldChange, padj= res$padj,genes=gsub("\\..*", "", rownames(res)))
plotdf<- plotdf %>%
  mutate(gene_type = case_when(log2FoldChange > 0 & pvalues < 0.05 ~ "upregulated",
                               log2FoldChange < 0 & pvalues < 0.05 ~ "downregulated",
                               TRUE ~ "not significant"))   
print(table(plotdf$gene_type))

dif.genes<- plotdf[which(plotdf$gene_type !="not significant"), ]


library(biomaRt)
mart <- useMart(dataset="hsapiens_gene_ensembl","ensembl")
G_list <- getBM(filters= "ensembl_gene_id", 
                attributes=c( 'ensembl_gene_id',
                              'hgnc_symbol', 'gene_biotype'),
                values= dif.genes$genes,
                mart= mart)

dif.genes<- left_join(dif.genes, G_list, by= c( "genes" = "ensembl_gene_id" ))


write.csv(dif.genes, file ="diff_genes_maleRUSH.csv", row.names = F)


