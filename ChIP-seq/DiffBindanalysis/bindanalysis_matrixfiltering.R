#This R files was used to create the matrices that are used for differential analysis inbetween the three groups young, old and AD

#read in file that got created by unionbedg
peak_df <- read.table("C:/Users/georg/Documents/Uni/Semester_2/Data Science in Life Sciences/Main_project/iCluster_prep/second_Laptop/diffBind_Rosa/H3K122ac_mat.txt")


#calculate length of the reads
length <- peak_df$V3 - peak_df$V2

#only plot for up to 2000 - to get an overview
hist(length[length<2000])

#longest read
max(length)


#check that in at least 23 of the 27 patients values are present so there is not too much missing values
sum(rowSums(peak_df[,4:30] != 0)>23)

#save the selected rows from previous analysis as new matrix
mat_23 <- peak_df[which(rowSums(peak_df[,4:30] != 0)>23),]


## functional analysis of the matrix - required packages
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
library(clusterProfiler)
library(ChIPpeakAnno)

#add index to GRanges object since somehow chromosome can not be accessed. -> to backtrack which region was selected as next to the TSS sites
mat_23$V31 <- 1:length(mat_23$V1)

#create GRanges object for finding out in next steps how far positions are from Transcription starting sites
gr <- toGRanges(mat_23, format="other",  colNames=c("space", "start", "end", "s1", "s2", "s3", "s4", "s5", "s6", "s7", "s8", "s9", "s10", "s11", "s12", "s13", "s14", "s15", "s16", "s17", "s18", "s19", "s20", "s21", "s22", "s23", "s24", "s25", "s26", "s27", "pos"))

promoter <- getPromoters(TxDb=txdb, upstream=3000, downstream=3000)

#create tagMatrix showing the distribution of the reads regarding distance to TSS
tagMatrix <- getTagMatrix(gr, windows=promoter)

tagHeatmap(tagMatrix, xlim=c(-3000, 3000), color="red")
summary(tagMatrix)

library(EnsDb.Hsapiens.v75) ##(hg19)
## create annotation file from EnsDb or TxDb
annoData <- toGRanges(EnsDb.Hsapiens.v75, feature="gene")


#try to get distance to next tss for each range
binOverFeature(gr, annotationData=annoData,
               radius=5000, nbins=20, FUN=mean, errFun=0,
               xlab="distance from TSS (bp)", ylab="count", 
               main="Distribution of aggregated peak numbers around TSS")

#only provides information about how near a range is to the next range
#tss <- distanceToNearest(gr)


overlaps.anno <- annotatePeakInBatch(gr, 
                                     AnnotationData=annoData, 
                                     output="nearestBiDirectionalPromoters",
                                     bindingRegion=c(-5000, 5000))

#show distribution of ranges in genome
pie1(table(overlaps.anno$insideFeature))

#show distribution of ranges in genome
out <- genomicElementDistribution(overlaps.anno, 
                                  TxDb = TxDb.Hsapiens.UCSC.hg19.knownGene,
                                  promoterRegion=c(upstream=5000, downstream=5000),
                                  geneDownstream=c(upstream=5000, downstream=5000),
                                  promoterLevel=list(
                                    # from 5' -> 3', fixed precedence 3' -> 5'
                                    breaks = c(-2000, -1000, -500, 0, 500),
                                    labels = c("upstream 1-2Kb", "upstream 0.5-1Kb", 
                                               "upstream <500b", "TSS - 500b"),
                                    colors = c("#FFE5CC", "#FFCA99", 
                                               "#FFAD65", "#FF8E32")))



#filter out the positions with TSS based on specified distances from the original matrix to obtain the exact positions in the genome
mat_TSS_filtered <- mat_23[overlaps.anno@elementMetadata@listData$pos,]

mat_TSS_filtered$TSS_dist <- overlaps.anno@elementMetadata$distance
mat_TSS_filtered$TSS_type <- overlaps.anno@elementMetadata$insideFeature
mat_TSS_filtered$gene_name <- overlaps.anno@elementMetadata$gene_name




count = 1
size = 1
resvec = c()
while(count < dim(mat_TSS_filtered)[1]){
  #check for correlation in columns next to each other
  while(ncol(intersect(mat_TSS_filtered[count,], mat_TSS_filtered[count+size,]))>21){
    size=size+1
  }
  #find column with least zeros (at the moment only first column taken)
  
  resvec <- append(resvec,count)
  
  count=count+size
  size=1
}  

#write created table into file, which then can be used for the actual binding analysis
write.table(mat_TSS_filtered[resvec,], file="C:/Users/georg/Documents/Uni/Semester_2/Data Science in Life Sciences/Main_project/iCluster_prep/second_Laptop/iCluster_mat_Rosa_H3K122ac.txt", row.names = FALSE, col.names = FALSE, sep="\t")










