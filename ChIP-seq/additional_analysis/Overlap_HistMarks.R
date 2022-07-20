# test Venn diagram for Young - AD 
library(ChIPpeakAnno)

# The MIT files are: (example for Y): all patient files of histone mark H3K9ac from Young patients are concatenated, sorted and merged = MIT file
file1 = "MIT.Y_H3K9ac_merged.bed"
Young <- toGRanges(file1, format="BED", header=FALSE)

file2 = "MIT.AD_H3K9ac_merged.bed"
AD <- toGRanges(file2, format="BED", header=FALSE)

file3 = "MIT.O_H3K9ac_merged.bed"
Old <- toGRanges(file3, format="BED", header=FALSE)

#find overlaps
ol3 <- findOverlapsOfPeaks(Young, AD, Old)


#create the Venn Diagram
V3 <- makeVennDiagram(ol3, fill=c("#CC79A7", "#56B4E9", "#F0E442"), col=c("#D55E00", "#0072B2", "#E69F00"), cat.col=c("#D55E00", "#0072B2", "#E69F00")) 











