#getwd()
#setwd("Dokumente/DataScience/MIT_files")

d1 <- read.table('MTL.H3K27ac.bed',header=F,sep="\t")
d1$peakwidth <- (d1[,3] - d1[,2])
d1$group <- "MTL.H3K27ac"
d1 <- subset(d1, select=c("peakwidth","group"))

d2 <- read.table('MTL.H3K9ac.bed',header=F,sep="\t")
d2$peakwidth <- (d2[,3] - d2[,2])
d2$group <- "MTL.H3K9ac"
d2 <- subset(d2, select=c("peakwidth","group"))

d3 <- read.table('MTL.H3K122ac.bed',header=F,sep="\t")
d3$peakwidth <- (d3[,3] - d3[,2])
d3$group <- "MTL.H3K122ac"
d3 <- subset(d3, select=c("peakwidth","group"))

d4 <- read.table('MTL.H3K4me1.bed',header=F,sep="\t")
d4$peakwidth <- (d4[,3] - d4[,2])
d4$group <- "MTL.H3K4me1"
d4 <- subset(d4, select=c("peakwidth","group"))



d <- rbind(d1,d2,d3,d4)
d$group <- factor(d$group,c("MTL.H3K27ac","MTL.H3K9ac","MTL.H3K122ac","MTL.H3K4me1"))


pdf('../peak_size.pdf',height=4,width=6)
ggplot(d,aes(log10(peakwidth),col=group)) + 
  #geom_density() + 
  geom_freqpoly() +
  theme_bw(base_size=16) +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  geom_hline(yintercept = 0)
dev.off()
