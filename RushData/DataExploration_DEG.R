## Differential Gene Expression Analysis For Rush 
library(tximport)
library(dplyr)
library(ggplot2)
library(DESeq2)
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
  
  overview <- mutate(overview, age.binary = age < 90 )
 

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

# create col/sample  data : condition sex, age of samples

condition<- factor(c(rep("C", length(files.CF) + length(files.CM) ),
                    rep("AD",length(files.ADF)+  length(files.ADM))
  
))

sex<- factor(c(rep("F", length(files.CF)),
               rep("M", length(files.CM)), 
               rep("F",length(files.ADF)),
               rep("M", length(files.ADM))
))

age.binary<- (c( overview.CF$age.binary,  overview.CM$age.binary, overview.ADF$age.binary, overview.ADM$age.binary ))
age<- (c( overview.CF$age,  overview.CM$age, overview.ADF$age, overview.ADM$age ))


sampletable <- data.frame(condition= condition, sex=sex, age.binary=age.binary, age=age)


# create dds object (summarizedexperiment obj) : design with binary age:  age >=90  vs. age < 90

dds.all <- DESeqDataSetFromTximport(txi = txi.rsem ,colData = sampletable ,design = ~ age.binary+ sex + condition)
dds <- DESeqDataSetFromTximport(txi = txi.rsem ,colData = sampletable ,design = ~ age.binary + condition)
#pre-filtering 
# ERCC counts, no genes
dds<- dds[1:(nrow(dds)-90),]
dds.all<- dds.all[1:(nrow(dds.all)-90),]
# low gene expression measurements across all samples
dds <- dds[ rowSums(counts(dds)) > 1, ]
dds.all <- dds.all[ rowSums(counts(dds.all)) > 1, ]

dds.female <- dds[, which(sex=="F")]
dds.male<- dds[, which(sex=="M")] 


# data exploration

# age distribution 
library(viridis)
library(forcats)
library(hrbrthemes)
library(ggridges)
ridges<- ggplot(sampletable, 
       aes(y =interaction(condition, sex), x= age, fill= interaction(condition, sex) )) +
  geom_density_ridges(alpha=0.6, stat="binline", bins=20) +
  theme_ridges() +
  theme(
    legend.position="none",
    panel.spacing = unit(0.1, "lines"),
    strip.text.x = element_text(size = 8)
  ) +
  xlab("") +
  ylab("condition and sex")

violins<- ggplot( sampletable, aes(x=interaction(condition, sex), y=age, fill=interaction( condition, sex))) +
  geom_violin(width=1.4) +
  geom_boxplot(width=0.1, color="black", alpha=0.2) +
  scale_fill_viridis(discrete = TRUE) +
  theme_ipsum() +
  theme(
    legend.position="none",
    plot.title = element_text(size=16),
    axis.title.y = element_text(hjust=0.55, vjust= 3 ,size= 13 ),
    axis.title.x = element_text(hjust=0.55 , vjust= -1.4 ,size= 13 )
    )+
  xlab("groups: condition + sex")





# variance stabilizing transformation. controlling great variance across the mean ( negative binomial data)
vsd<- vst(dds.all, blind=FALSE)
# or  regularized-logarithm transformation, takes more time
#rld <- rlog(dds, blind = FALSE)




dds.all <- estimateSizeFactors(dds.all)


create_log2_comparsion <- function(rows){
  
  counts <- counts(dds.all, normalized=T)
  colnames<- colnames(counts)[rows]
  df <- bind_rows(
    as_data_frame(log2(counts[, rows]+1)) %>%
      mutate(transformation = "log2(x + 1)"),
    as_data_frame(assay(vsd)[, rows]) %>% mutate(transformation = "vst"),
  #  as_data_frame(assay(rld)[, 1:2]) %>% mutate(transformation = "rlog")
  )
  
  
  colnames(df)[1:2] <- c("x", "y")  
  
  plot<- ggplot(df, aes(x = x, y = y))+ geom_hex(bins = 80) +
  coord_fixed() + facet_grid( . ~ transformation)  +
    xlab("") +
    ylab("")
  
  return(plot)
}


## random samples
plot_cf<- create_log2_comparsion(c(1,2))
plot_cm<- create_log2_comparsion(c(27,28))
plot_control<- create_log2_comparsion(c(1, 27))
plot_af <- create_log2_comparsion(c(42,43))
plot_am <- create_log2_comparsion(c(71,72))
plot_AD<- create_log2_comparsion(c(42,71))

## samples with highest/lowest correlation 
cor.mat<- counts(dds.all , normalized=T) %>% cor( method="spearman")
min.cor <- which(cor.mat ==min(cor.mat), arr.ind = T)
min.val <- min(cor.mat)
cor.mat[which(cor.mat ==max(cor.mat), arr.ind = T)]=0
max.cor <-which(cor.mat ==max(cor.mat), arr.ind = T)
max.val<- max(cor.mat)


# highest correlation between two samples of group: control + male
# lowest correaltion between sample of AD male and Control Female 
plot.max.cor <- create_log2_comparsion(max.cor[1,]) + ylab(paste("highest correlation: ", round(max.val, digit=3)))
plot.min.cor <- create_log2_comparsion(min.cor[1,]) + ylab(paste("lowest correlation: ", round(min.val, digit=3)))


# sample sample distance 

#PCA

pcaData <- plotPCA(vsd, intgroup= c("condition", "age.binary", "sex"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
pcaplot<- ggplot(pcaData, aes(PC1, PC2, shape=sex, color=condition)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) 


# heatmap 

library("pheatmap")
library("RColorBrewer")
sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd$condition, vsd$sex, vsd$age, sep="-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pdf("heatmap.pdf", height=10)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)
dev.off()



### differential gene expression analysis
analyzed.dds.all <- DESeq(dds.all)
analyzed.dds.male <- DESeq(dds.male)
analyzed.dds.female <- DESeq(dds.female)

res.dds.all <- results(analyzed.dds.all ,alpha=0.05)
res.dds.male <- results(analyzed.dds.male ,alpha=0.05)
res.dds.female <- results(analyzed.dds.female ,alpha=0.05)


res<- res.dds.all

resOrdered <- res[order(res$pvalue),]

summary(res)

pvalues <- p.adjust(res$pvalue, method = "BH")
sum(pvalues<0.05,na.rm = T)

# only 7  diff. expr. gene detected !!!

# MA plot
plotMA(res, ylim=c(-4,4))

resNorm <- lfcShrink(analyzed.dds.all, coef=2, type="normal")
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


write.csv(dif.genes, file ="diff_genes_allRUSH.csv", row.names = F)



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



