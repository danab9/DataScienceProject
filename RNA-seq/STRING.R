#STRING ANALYSIS
library(stringr)
library(ggplot2)

rnaseqdata <- read.table("RNAseq/normalizedRNA.txt")

# import genes from GO enrichment "Regulation of transcription"
ownres <- unlist(read.table( 'go_genes.txt'))
ownres <- ownres[!ownres%in%c('LINC-PINT','SGF29','TFDP2','ZNF137P')]


# genexpression for color
upregulated <- read.table("upregulated.txt",header=TRUE)
stringup <- upregulated[upregulated$genes%in%ownres,]              
stringup <- stringup[order(stringup$genes),]
round(stringup$log2FoldChange*100)

#expression Values for node size
AD <-  dplyr::select(rnaseqdata, starts_with("AD"))
stringADnorm <- AD[order(rownames(AD)),]
rowMeans(stringADnorm)*10

edges <- read.csv("STRING_network_default_edge.csv", header=TRUE, stringsAsFactors=FALSE)
names <- unlist(strsplit(edges$name,"[ (pp) ]"))
names <- matrix(names[names!=""],ncol=2,byrow=TRUE)
order <- sort(table(unlist(names)),decreasing = TRUE)
data <- data.frame(order)

# Basic barplot
png(file = "String_interactions.png", height = 400,width=600)
p<-ggplot(data=data[1:13,], aes(x=Var1, y=Freq)) +
  geom_bar(stat="identity", fill="grey30")+
  geom_text(aes(label=Freq), vjust=-0.3, size=4.5)+
  theme_minimal()
p + labs(title="", 
     x="", y = "# String interactions",size=5) +
  theme(axis.text=element_text(size=11),
        axis.title=element_text(size=13),
        axis.text.x = element_text(angle = 45, vjust = 0.7, hjust=0.3,size=12))
dev.off()
