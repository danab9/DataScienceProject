#getwd()
#setwd("Dokumente/DataScience/MIT_files/MTL_files(ChIPS)")

d1 <- read.table('multiMTL.bed',header=F,sep="\t") #read in the MTL files (=MIT files of three different patient groups with same histone mark)
d1$peakwidth <- (d1[,3] - d1[,2]) #check for length of peaks
d1$group <- "multiMTL"
d1 <- subset(d1, select=c("peakwidth","group"))

d2 <- read.table('MTL.H3K27ac_merged.bed',header=F,sep="\t")
d2$peakwidth <- (d2[,3] - d2[,2])
d2$group <- "MTL.H3K27ac"
d2 <- subset(d2, select=c("peakwidth","group"))

d3 <- read.table('MTL.H3K9ac_merged.bed',header=F,sep="\t")
d3$peakwidth <- (d3[,3] - d3[,2])
d3$group <- "MTL.H3K9ac"
d3 <- subset(d3, select=c("peakwidth","group"))

d4 <- read.table('MTL.H3K122ac_merged.bed',header=F,sep="\t")
d4$peakwidth <- (d4[,3] - d4[,2])
d4$group <- "MTL.H3K122ac"
d4 <- subset(d4, select=c("peakwidth","group"))




d <- rbind(d1,d2,d3,d4)
d$group <- factor(d$group,c("multiMTL", "MTL.H3K27ac","MTL.H3K9ac","MTL.H3K122ac"))


# create and save the peak sizes of the three different acetylation MTLs and multiMTL file
pdf('../peak_size.pdf',height=4,width=6)
ggplot(d,aes(log10(peakwidth),col=group)) + 
  #geom_density() + 
  geom_freqpoly() +
  theme_bw(base_size=16) +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  geom_hline(yintercept = 0)
dev.off()
