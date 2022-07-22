### For creation of Matrix, which is required for iClusterBayes --- !!! add code before (from second laptop)

mat <- read.table("iCluster_matrix.bed")

#check that 
sum(rowSums(mat != 0)>=27)
mat_27 <- mat[which(rowSums(mat != 0)>=27),]
range_length2 = mat_27$V3 - mat_27$V2
hist(range_length2[range_length2 < 15])
sum(range_length2>15)

mat_new <- mat_27[which(range_length2>15),]


# check that ranges are not too big

range_length = mat_new$V3 - mat_new$V2
hist(range_length[range_length < 15])
sum(range_length>15)
sum(range_length<500)
mat_500 <- mat_new[which(range_length<500),]



## functional analysis of the matrix
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
library(clusterProfiler)
library(ChIPpeakAnno)

#add index to GRanges object since somehow chromosome can not be accessed. -> to backtrack which region was selected as next to the TSS sites
mat_500$V31 <- 1:length(mat_500$V1)



gr <- toGRanges(mat_500, format="other", space=mat_500$V1, colNames=c("space", "start", "end", "s1", "s2", "s3", "s4", "s5", "s6", "s7", "s8", "s9", "s10", "s11", "s12", "s13", "s14", "s15", "s16", "s17", "s18", "s19", "s20", "s21", "s22", "s23", "s24", "s25", "s26", "s27", "pos"))
#covplot(gr, weightCol = "s1")

promoter <- getPromoters(TxDb=txdb, upstream=3000, downstream=3000)
tagMatrix <- getTagMatrix(gr, windows=promoter)

#create matrix for reads with max length of 500
tagHeatmap(tagMatrix, xlim=c(-3000, 3000), color="red")
summary(tagMatrix)


library(EnsDb.Hsapiens.v75) ##(hg19)
## create annotation file from EnsDb or TxDb
annoData <- toGRanges(EnsDb.Hsapiens.v75, feature="gene")


#try to get distance to next tss for each range
binOverFeature(gr, annotationData=annoData,
               radius=5000, nbins=20, FUN=length, errFun=0,
               xlab="distance from TSS (bp)", ylab="count", 
               main="Distribution of aggregated peak numbers around TSS")

#only provides information about how near a range is to the next range
#tss <- distanceToNearest(gr)


overlaps.anno <- annotatePeakInBatch(gr, 
                                     AnnotationData=annoData, 
                                     output="nearestBiDirectionalPromoters",
                                     bindingRegion=c(-2000, 500))

#show distribution of ranges in genome
pie1(table(overlaps.anno$insideFeature))

#show distribution of ranges in genome
out <- genomicElementDistribution(overlaps.anno, 
                                  TxDb = TxDb.Hsapiens.UCSC.hg19.knownGene,
                                  promoterRegion=c(upstream=2000, downstream=500),
                                  geneDownstream=c(upstream=0, downstream=2000),
                                  promoterLevel=list(
                                    # from 5' -> 3', fixed precedence 3' -> 5'
                                    breaks = c(-2000, -1000, -500, 0, 500),
                                    labels = c("upstream 1-2Kb", "upstream 0.5-1Kb", 
                                               "upstream <500b", "TSS - 500b"),
                                    colors = c("#FFE5CC", "#FFCA99", 
                                               "#FFAD65", "#FF8E32")))



#filter out the positions with TSS based on specified distances from the original matrix to obtain the exact positions in the genome
mat_TSS_filtered <- mat_500[overlaps.anno@elementMetadata@listData$pos,]

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
  #print(size)
  #find column with least zeros (at the moment only first column taken)
  
  #print(size)
  #print(count)
  resvec <- append(resvec,count)
  
  count=count+size
  size=1
  #print("check")
}  
df_undoubled <- mat_TSS_filtered[resvec,]

write.table(df_undoubled, file="iCluster_mat_12k.txt", row.names = FALSE, col.names = FALSE, sep="\t")

