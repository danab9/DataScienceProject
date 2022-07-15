## differential enrichdfent analysis
## codfparing nordfalized signal between AD and old for each peak 
## dfethod like in paper: Wilc. Rank test , two-sided, pval =0.05 

df <- read.table("~/Downloads/iCluster_mat_12k_Rosa.txt")
AD_index <- c(22:30)
old_index <- c(12:21)
young_index <- c(4:11)


pvals <- c()
foldchanges<- c()
for (i in 1: nrow(df)){
  grp2 <- as.numeric(df[i, AD_index])
  grp1<- as.numeric( df[i, old_index])
  fc<- mean(grp2) -mean(grp1)
  log2FC<- sign(fc)* log2(abs(fc))
  p<- wilcox.test(grp1, grp2)$p.value
  pvals[i] <- p
  foldchanges[i] <- log2FC
  
}

colnames(df)[32] <- "distance.to.TSS"
colnames(df)[34] <- "gene.name"
colnames(df)[33] <- "region"


df$pval <- as.numeric(pvals)
df$log2FC <- as.numeric(foldchanges)

adj.pval<- p.adjust(df$pval, method="BH")
df$adj.pval <- adj.pval


df<- df %>%
  mutate(status = case_when(log2FC > 0 & pval < 0.05 ~ "gain",
                               log2FC < 0 & pval < 0.05 ~ "loss",
 
                              TRUE ~ "not significant"))   

df<- df %>%
  mutate(status2 = case_when( pval < 0.05 ~ "significant",
              
                            TRUE ~ "not significant"))   



df<- df %>%
  mutate(TSS = case_when( distance.to.TSS <5 ~ "<5",
                          distance.to.TSS > 5  & distance.to.TSS <= 25 ~ "5-25",
                          distance.to.TSS > 25  & distance.to.TSS <= 50 ~ "25-50",
                          distance.to.TSS > 50  & distance.to.TSS <= 100 ~ "50-100",
                          distance.to.TSS > 100 & distance.to.TSS <= 300 ~ "100-300",
                          distance.to.TSS > 300 ~ ">300",
                          
                          ) )   




# volcanoplot
ggplot(df, aes(y= - log10(pval), x=log2FC, color=status2)) + geom_point() + xlim(-10,10)

# histogram
ggplot(filter(df), aes(x=log2FC, color=status2)) + geom_density()
  geom_histogram()


g<-filter(df, status!="not significant" ) %>% group_by(TSS) %>% summarise(count=sum(status=="gain")) %>% cbind(status=rep("gain", 6))
l<- filter(df, status!="not significant" ) %>% group_by(TSS) %>% summarise(count=sum(status=="loss")) %>% cbind(status=rep("loss", 6))
plot.df <- rbind(g,l)
plot.df$TSS <- factor(plot.df$TSS, levels = c("<5", "5-25", "25-50", "50-100", "100-300", ">300"))   
    
ggplot(plot.df, aes(x = TSS, y=count ,fill = status)) + 
  geom_bar(stat="identity", color = "black") +
  scale_fill_manual(values = c( "#9E9AC8", "#6A51A3")) +
  guides(fill = guide_legend(title = "")) +
  ylab("") 

## ANOVA test

## checking normality 
library(car)

shapiro.res<- c()
levene.res<- c()
anova.res <- c()
for (i in 1:nrow(df)){
  grp1 <- data.frame(count=as.numeric(df[i, AD_index]), group=rep("AD", length(AD_index)))
  grp2 <- data.frame(count=as.numeric(df[i, old_index]), group=rep("old", length(old_index)))
  grp3 <- data.frame(count=as.numeric(df[i, young_index]), group=rep("young", length(young_index)))

  helper<- rbind(grp1,grp2,grp3)

  res_aov <- aov(count ~ group, data=helper)
  
  anova.res[i] <- unlist(summary(res_aov))[9]
  #shapiro.res[i]<- shapiro.test(res_aov$residuals)$p.value
  #levene.res[i] <- leveneTest(count ~group, data=helper)$Pr[1]
}

length(which(anova.res < 0.05))
#adj.anova<- p.adjust(anova.res, method="BH")
# qqplot<- qqPlot(res_aov$residuals, id = FALSE )

df <- cbind(df, anova.pval = anova.res)

# intersection of anova and wilcoxon test

disease.specific<- filter(df, pval <0.05 & anova.pval <0.05)

##boxplot 

gain.dis.spec <- filter(disease.specific, status=="gain")
loss.dis.spec <- filter(disease.specific, status=="loss")

loss.all.young <-unlist(loss.dis.spec[, young_index])
loss.all.old <- unlist(loss.dis.spec[, old_index])
loss.all.AD <- unlist(loss.dis.spec[, AD_index])
loss.labels<- c( rep("young", length(loss.all.young)), 
            rep("old", length(loss.all.old)),
            rep("AD", length(loss.all.AD))
            
)



all.young <-unlist(gain.dis.spec[, young_index])
all.old <- unlist(gain.dis.spec[, old_index])
all.AD <- unlist(gain.dis.spec[, AD_index])
labels<- c( rep("young", length(all.young)), 
            rep("old", length(all.old)),
            rep("AD", length(all.AD))
  
)

gain.length <- length(all.young)+ length(all.old)+ length(all.AD)
loss.length <- length(loss.all.young)+ length(loss.all.old)+ length(loss.all.AD)


boxplot.df <- data.frame(counts= c(all.young, all.old, all.AD, loss.all.young, loss.all.old, loss.all.AD), 
                         group=c(labels, loss.labels),
                         status= c(rep("gain", gain.length ),rep("loss", loss.length ) ))


ggplot(boxplot.df, aes(x=interaction(group,status), y=counts, color=group)) +
  geom_boxplot()



write.csv(disease.specific, file= "~/disease_spec_enriched_peaks.csv", row.names = F)


# comparing it to DEG genes
down_DEG <- read.table("/home/rosa/DataScienceProject/RNA-seq/UpDownreg_genes/downregulationgenes.txt")
up_DEG <- read.table("/home/rosa/DataScienceProject/RNA-seq/UpDownreg_genes/upregulatedgenes.txt")

down.inter.genes<- intersect(down_DEG$V1, disease.specific$gene.name)
up.inter.genes<- intersect(up_DEG$V1, disease.specific$gene.name)

gene.label <- c(rep("down", length(down.inter.genes)), rep("up", length(up.inter.genes)))

down.status<- filter(disease.specific, gene.name %in% down.inter.genes)$status 
up.status<- filter(disease.specific, gene.name %in% up.inter.genes)$status 
status.label <- c(down.status, up.status)

intersection.df <- data.frame(gene.name = c(down.inter.genes, up.inter.genes),
                              DEG.status= gene.label,
                              chip.status = status.label
                                
                                )

write.csv(intersection.df, file="~/intersection_rnaseq_chipseq.csv", row.names = F)
