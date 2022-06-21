library(DESeq2)
# http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html
diffexpression <- function(group1,group2,group3=NA,alpha=0.05,tidy=TRUE,result=TRUE){
  if(all(is.na(group3))){
    condition <- factor(c(rep("G1",ncol(group1)),rep("G2",ncol(group2))))
    dds <- DESeqDataSetFromMatrix(countData = data.frame(group1,group2),
                                  DataFrame(condition), ~ condition)
  }
  else{
    condition <- factor(c(rep("G1",ncol(group1)),rep("G2",ncol(group2)),rep("G3",ncol(group3))))
    dds <- DESeqDataSetFromMatrix(countData = data.frame(group1,group2,group3),
                                  DataFrame(condition), ~ condition)
  }
  dds <- DESeq(dds)
  if(!result){
    return(dds) # to obtain also the genes
  }
  return(results(dds,alpha=alpha,tidy=tidy))
}
