library(DESeq2)
library(dplyr)

# import data
rnaseq_data <- read.csv("/home/rosa/Dokumente/GSE153873_summary_count.star.txt", sep = "\t")
row.names(rnaseq_data)<- rnaseq_data[,1]

# make subset of AD, young, and old
AD <- select(rnaseq_data, ends_with("AD"))
young<- select(rnaseq_data, ends_with("Young"))
old<- select(rnaseq_data, ends_with("Old"))

apply(old, 2, function(x) any(is.na(x)))

# function to make DEG
make_DGE<- function(group1, group2, test){
  cond <- data.frame(sample=c(colnames(group1), colnames(group2)), condition= c(rep(0,ncol(group1)),rep(1,ncol(group2))))
  
  group1 <- cbind(refGene=rownames(group1), group1)
  group2 <- cbind(refGene=rownames(group2), group2)

  count_data<- select(inner_join(group1, group2, by="refGene"), -refGene)
  #add pseudo counts count_data <- count_data +1
  print(apply(count_data, 2, function(x) any(is.na(x))))
  
  if (test == "DESeq2"){
    dds<- DESeqDataSetFromMatrix(countData= count_data ,colData=cond, design=~ condition)
    dds<- DESeq(dds)
    res<- results(dds)
    return(res)
  }
 
}

deseq_ad_old <- make_DGE(AD, old , "DESeq2")



library(ggplot2)



make_volcanoplot <- function(results, path, name ){
  
  plotdf<- data_frame(log2FoldChange=results$log2FoldChange, padj= results$padj)
  
  plotdf<- plotdf %>%
    mutate(gene_type = case_when(log2FoldChange >= log2(2) & padj <= 0.05 ~ "upregulated",
                                 log2FoldChange <= log2(0.5) & padj <= 0.05 ~ "downregulated",
                                 TRUE ~ "not significant"))   
  
  cols <- c("upregulated" = "#ffad73", "downregulated" = "#26b3ff", "not significant" = "grey") 
  sizes <- c("upregulated" = 1.5, "downregulated" = 1.5, "not significant" = 1) 
  alphas <- c("upregulated" = 1, "downregulated" = 1, "not significant" = 0.5)
  
  
  volcanplot<- ggplot(plotdf, aes(x = log2FoldChange,
                     y =  -(log10(padj)),
                     fill = gene_type,    
                     size = gene_type,
                     alpha = gene_type
  )) + 
    geom_point(shape = 21, # Specify shape and colour as fixed local parameters    
               colour = "black") + 
    geom_hline(yintercept = -log10(0.05),
               linetype = "dashed") + 
    geom_vline(xintercept = c(log2(0.5), log2(2)),
               linetype = "dashed") +
    scale_fill_manual(values = cols) + # Modify point colour
    scale_size_manual(values = sizes) + # Modify point size
    scale_alpha_manual(values = alphas)  # Modify point transparency
   
  ggsave(filename= name, volcanplot, path = path)
  
  print(table(plotdf$gene_type))
}



path = "/home/rosa/DataScienceProject"


make_volcanoplot(deseq_ad_old, path, "volcanplot.png")

 


