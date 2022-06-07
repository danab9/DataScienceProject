library(DESeq2)
library(dplyr)

# import data
rnaseq_data <- read.csv("/home/rosa/Dokumente/GSE153873_summary_count.star.txt", sep = "\t")
row.names(rnaseq_data)<- rnaseq_data[,1]
rnaseq_data<- select(rnaseq_data, -refGene)

# make subset of AD, young, and old
AD <- select(rnaseq_data, ends_with("AD"))
young<- select(rnaseq_data, ends_with("Young"))
old<- select(rnaseq_data, ends_with("Old"))


# function to make DEG
make_DGE<- function(group1, group2, test){
  cond <- data.frame(sample=c(colnames(group1), colnames(group2)), condition= c(rep(0,ncol(group1)),rep(1,ncol(group2))))
  count_data<- inner_join(group1, group2)
  if (test == "DESeq2"){
    dds<- DESeqDataSetFromMatrix(countData= count_data ,colData=cond, design=~ condition)
    dds<- DESeq(dds)
    res<- results(dds)
    return(res)
  }
 
}

deseq_ad_old <- make_DGE(AD, old , "DESeq2")

