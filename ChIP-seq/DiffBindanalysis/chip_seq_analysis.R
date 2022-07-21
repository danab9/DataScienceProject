## differential binding analysis
## Since the diffBind R-package didnt work due to missing qualitative bam files, we followed the methods of Nativio et al. 
## we are comparing the signal values between the groups for each peak/region 
## those signals are already normalized by RPKM ( reads per kilobase millions)
## we used a Wilcoxon Rank sum test for the pairwise comparison and for a further validation ANOVA for multiple comparison 

## the following script is repeated for each histone modificaiton data set

## reading in 
df <- read.table("~/DataScienceProject/data/mat_H3K9ac.txt")
AD_index <- c(22:30)
old_index <- c(12:21)
young_index <- c(4:11)

colnames(df)[32] <- "distance.to.TSS"
colnames(df)[34] <- "gene.name"
colnames(df)[33] <- "region"

# function for the wilcoxon test 
# no p-value adjustment ! 
execute_wilcox_test <- function(df, name, index1,index2){
  name.pv <- paste(name, "pval", sep=".")
  name.fc <- paste(name, "log2FC", sep=".")
  pvals <- c()
  foldchanges<- c()
  for (i in 1: nrow(df)){
    grp2 <- as.numeric(df[i, index2])
    grp1<- as.numeric( df[i, index1])
    fc<- (mean(grp2)/mean(grp1))
    log2FC<- log2(fc)
    p<- wilcox.test(grp1, grp2)$p.value
    pvals[i] <- p
    foldchanges[i] <- log2FC
    
  }
  df[, name.pv] <- pvals
  df[, name.fc ]<- foldchanges
  
  
  if (name =="ADO"){
    df<- df %>%
      mutate(ADO.status = case_when(ADO.log2FC > 0 & ADO.pval < 0.05 ~ "gain",
                                     ADO.log2FC < 0 & ADO.pval < 0.05 ~ "loss",
                                     TRUE ~ "not significant"))  
    print(table(df$ADO.status))
    
  }
  else{
    df<- df %>%
      mutate(OY.status = case_when(OY.log2FC > 0 & OY.pval < 0.05 ~ "gain",
                                     OY.log2FC < 0 & OY.pval < 0.05 ~ "loss",
                                     
                                     TRUE ~ "not significant"))   
    print(table(df$OY.status))
    
  }
 
  return(df)
}


df <- execute_wilcox_test(df, "ADO", old_index, AD_index)
df <- execute_wilcox_test(df, "OY", young_index, old_index)




## plotting results


# 1. volcano plot , visualizing the distribution of the log2FC vs p-values, 
# for H3K27ac and H3K9ac more gains than losses, H3K122ac the other way around
create_volcanoplot <- function(df, name){
  name.pval <- paste(name, "pval", sep=".")
  name.fc <- paste(name, "log2FC", sep=".")
  data <-df[, c(name.pval, name.fc)]
  colnames(data)<- c("pval", "log2FC")  
  data <- data %>%
      mutate(status2 = case_when( pval < 0.05 ~ "significant",
                                  
                                  TRUE ~ "not significant"))   
    
  plot <- ggplot(data, aes(y= - log10(pval), x=log2FC, color=status2)) + geom_point()+ labs(colour=NULL)
  return (plot)
  
}

create_volcanoplot(df, "ADO")

#2. distribution of log2FC divided into not significant regions and significant ones
# to have better look on log2 FC wothout p-vals, non significant ones are quite normally distributed
create_FC_dist <- function(df, name, type="density"){
  name.pval <- paste(name, "pval", sep=".")
  name.fc <- paste(name, "log2FC", sep=".")
  data <-df[, c(name.pval, name.fc)]
  colnames(data)<- c("pval", "log2FC")  
  data <- data %>%
    mutate(status2 = case_when( pval < 0.05 ~ "significant",
                                
                                TRUE ~ "not significant"))   
  
  if (type == "density"){
    return( ggplot(data, aes(x=log2FC, color=status2)) + geom_density() + labs(color=NULL))
  }
  else {
    return(ggplot(data, aes(x=log2FC, color=status2)) +   geom_histogram() + labs(color=NULL))
  }
  
}

create_FC_dist(df, "ADO", "hist") + xlim(-5,5) + ylab("")


#3. Distribution of significant peaks according to their distances to the next TSS (Transcription Starting Site)
# most of the peaks are either very close to the TSS or far away, ths result differs from Nativio et al outcome

create_barchart_TSS <- function(df, name){
  name.status <- paste(name, "status", sep=".")
  data <-df[, c(name.status, "distance.to.TSS")]
  colnames(data)<- c("status", "distance.to.TSS")  
  data<- data %>%
    mutate(TSS = case_when( distance.to.TSS <=5 ~ "<5",
                            distance.to.TSS > 5  & distance.to.TSS <= 100 ~ "5-100",
                            distance.to.TSS > 100 & distance.to.TSS <= 500 ~ "100-500",
                            distance.to.TSS > 500 & distance.to.TSS <= 1000 ~ "500-1000",
                            distance.to.TSS > 1000 ~ ">1000"
                            
    ) )  
  data$TSS <- factor(data$TSS, levels = c("<5", "5-100", "100-500", "500-1000", ">1000"))   
  
  gains<-filter(data, status!="not significant" ) %>% group_by(TSS) %>% summarise(count=sum(status=="gain")) %>% cbind(status=rep("gain", 5))
  losses<- filter(data, status!="not significant" ) %>% group_by(TSS) %>% summarise(count=sum(status=="loss")) %>% cbind(status=rep("loss", 5))
  plot.df <- rbind(gains,losses)
  
  plot <- ggplot(plot.df, aes(x = TSS, y=count ,fill = status)) + 
    geom_bar(stat="identity", color = "black") +
    scale_fill_manual(values = c( "#9E9AC8", "#6A51A3")) +
    guides(fill = guide_legend(title = "")) +
    ylab("") 
  
  return(plot)
  
}

create_barchart_TSS(df, "ADO")


## ANOVA test
# to validate result of Wilcoxon test

execute.ANOVA.test <- function(df, idx1,idx2,idx3){
  ## checking assumption of normality and equal sample variances
  
  shapiro.res<- c()
  levene.res<- c()
  anova.res <- c()
  for (i in 1:nrow(df)){
    grp1 <- data.frame(count=as.numeric(df[i, idx1]), group=rep("AD", length(idx1)))
    grp2 <- data.frame(count=as.numeric(df[i, idx2]), group=rep("old", length(idx2)))
    grp3 <- data.frame(count=as.numeric(df[i, idx3]), group=rep("young", length(idx3)))
    
    together<- rbind(grp1,grp2,grp3)
    
    res_aov <- aov(count ~ group, data=together)
    
    anova.res[i] <- unlist(summary(res_aov))[9]
    shapiro.res[i]<- shapiro.test(res_aov$residuals)$p.value
    levene.res[i] <- leveneTest(count ~group, data=together)$Pr[1]
  }
  
  results <- data.frame( anova.pval = anova.res, shapiro.pval = shapiro.res, levene.pval = levene.res )
  
  return(results)
}


### plotting to check normality , alternative to the shapiro test
  #library(car)
  # qqplot<- qqPlot(res_aov$residuals, id = FALSE )

anova.df <- execute.ANOVA.test(df, young_index, old_index, AD_index)

# normality is not given for most cases! 
failed.normality.test <- mean(anova.df$shapiro.pval < 0.05)
failed.eq.sample.variance <- mean(anova.df$levene.pval < 0.05)


# intersection of anova and wilcoxon test
df <- cbind(df, anova.pval=anova.df$anova.pval)

disease.specific<- filter(df, ADO.pval <0.05 & anova.pval <0.05)
age.regulated <- filter(df,OY.pval <0.05 & ADO.pval >= 0.05 & anova.pval <0.05 )
age.dysregulated <- filter(df, (ADO.pval <0.05 | OY.pval <0.05) & anova.pval <0.05)

##boxplot
## have a look on the overall distribution of the signals for the significant peaks between the groups 
# as expected, AD has a lower median for losses and a higher median for gains compared to the old and young group 
create.boxplot <- function(data, name ){
  name.status <- paste(name, "status", sep=".")
  
  gain.df<- data[which(data[, name.status] == "gain"),]
  loss.df<- data[which(data[, name.status] == "loss"),]
  
  loss.all.young <-unlist(loss.df[, young_index])
  loss.all.old <- unlist(loss.df[, old_index])
  loss.all.AD <- unlist(loss.df[, AD_index])
  loss.labels<- c( rep("young", length(loss.all.young)), 
                   rep("old", length(loss.all.old)),
                   rep("AD", length(loss.all.AD))
                   
  )
  
  all.young <-unlist(gain.df[, young_index])
  all.old <- unlist(gain.df[, old_index])
  all.AD <- unlist(gain.df[, AD_index])
  labels<- c( rep("young", length(all.young)), 
              rep("old", length(all.old)),
              rep("AD", length(all.AD))
              
  )
  
  gain.length <- length(all.young)+ length(all.old)+ length(all.AD)
  loss.length <- length(loss.all.young)+ length(loss.all.old)+ length(loss.all.AD)
  
  boxplot.df <- data.frame(counts= c(all.young, all.old, all.AD, loss.all.young, loss.all.old, loss.all.AD), 
                           group=c(labels, loss.labels),
                           status= c(rep("gain", gain.length ),rep("loss", loss.length ) ))
  
  
  plot<- ggplot(boxplot.df, aes(x=group, y=counts, color=group)) +
    geom_boxplot()  + 
    theme( axis.text.x = element_blank())+   
    xlab("")+
    ylab("")+
    facet_grid(~status)
  
  return(plot)
  
  
}


create.boxplot(disease.specific, "ADO")

# saving all disease specific differentially binded regions, all three files are saved in the folder results
write.csv(disease.specific, file= "~/H3K9ac_disease_spec_enriched_peaks.csv", row.names = F)



#  CHip-Seq and RNA-seq

# comparing the peaks and their associated genes with the DEG genes of the RNA-seq 
# making venn diagrams to see the intersections 

# reading in the significant DEG
library(ggvenn)
down.path <- "~/DataScienceProject/RNA-seq/UpDownreg_genes/downregulationgenes.txt"
up.path <- "~/DataScienceProject/RNA-seq/UpDownreg_genes/upregulatedgenes.txt"
save.path <- "~/intersection_H3K9ac.csv"

# creating intersections, venn diagrams and saving everything in a table
# we expect intersections between up DEGs and gains, since the histone acetylation are enhancing transcription 
# and we expect intersections between down DEGs and losses, since with lossing acetylation the trnascription is not activated. 
# however, we oberserved also other intersections, which makes biologically not a lot of sense, and in genereal a rather smaller overlap between the two omic data sets. 
comapre_chip_rnaseq  <- function(data, up.path, down.path, save.path){
  
  down_DEG <- read.table(down.path)
  up_DEG <- read.table(up.path)
  
  down.inter.genes<- intersect(down_DEG$V1, disease.specific$gene.name)
  up.inter.genes<- intersect(up_DEG$V1, disease.specific$gene.name)
  
  gene.label <- c(rep("down", length(down.inter.genes)), rep("up", length(up.inter.genes)))
  
  down.status<- filter(filter(disease.specific, gene.name %in% down.inter.genes) , !duplicated(gene.name))$ADO.status 
  up.status<- filter(filter(disease.specific, gene.name %in% up.inter.genes), !duplicated(gene.name))$ADO.status 
  
  status.label <- c(down.status, up.status)
  
  gain.gene.names <- filter(filter(disease.specific, ADO.status=="gain") , !duplicated(gene.name))$gene.name
  loss.gene.names <- filter(filter(disease.specific, ADO.status=="loss") , !duplicated(gene.name))$gene.name
  
  
  
  intersection.df <- data.frame(gene.name = c(down.inter.genes, up.inter.genes),
                                DEG.status= gene.label,
                                chip.status = status.label
                                
  )
  
  write.csv(intersection.df, file=save.path, row.names = F)
  

  plot<- ggvenn(list(downDEG=down_DEG$V1 , upDEG=up_DEG$V1 ,loss=loss.gene.names,gain=gain.gene.names  ))
  
  return(plot)
}

comapre_chip_rnaseq(disease.specific, up.path, down.path , save.path)




