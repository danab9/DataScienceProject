# icluster
library(iClusterPlus)
library(dplyr)
# get k:
bayfit = tune.iClusterBayes(cpus=1, dt1=rnaseq_all, dt2=data.matrix(chipseq_12k),
                            type=c("poisson","poisson"), K=1:7)

allBIC = NULL 
devratio=NULL
nK = length(bayfit$fit)
for(i in 1:nK){
  allBIC = c(allBIC, bayfit$fit[[i]]$dev.ratio)
  devratio = c(devratio, bayfit$fit[[i]]$dev.ratio)
}
# save plots of Ks performance  # emphasized: best optional Ks to notice.
png("plots/iCluster_Ks.png")
par(mar=c(4.0,4.0, 0.5,0.5), mfrow=c(1,2))
plot(1:nK, allBIC,type="b", xlab="l", ylab="BIC", pch=c(1,1,19,19,19,1,1))
plot(1:nK, devratio, type="b", xlab="k",ylab="Deviance ratio", pch=c(1,1,19,19,19,1,1))
dev.off()

# run iCluster Bayes with chosen K (K=4)
iclust_res_all <- iClusterBayes(dt1=rnaseq_all,dt2=data.matrix(chipseq_12k), type = c("poisson", "poisson"),K=4)

# save iCluster figs as png
png("plots/iClusterBayes_K4.png")
plotHMBayes(fit=iclust_res_all, datasets=list(rnaseq_all,data.matrix(chipseq_12k)),type=c("poisson", "poisson"), sparse = c(1,1), scale = c(T,F))
dev.off()

#clusters <- iclust_res_all$clusters
#length(clusters)


# load metadata
metadata <- read.csv(file="data/metadata.tsv", sep='\t')  # To include 
row.names(metadata) <- paste(metadata$Study.Group,metadata$Sample.ID, sep='X')
metadata <- metadata[-which(!row.names(metadata) %in% row.names(rnaseq_all)),]
metadata <- metadata %>% select(-c(Sample.ID))
#dim(metadata)
metadata <- metadata[order(row.names(metadata)),]


# AD vs. control
library(mixOmics)
y <- setNames(metadata$Study.Group, row.names(metadata))
y[which(y!="AD")] <- "control"
y <- y[order(names(y))]
y <- y[which(names(y) %in% row.names(rnaseq_all))]
# make sure chipseq names are unique
chipnames <- colnames(chipseq_12k)
colnames(chipseq_12k) <- make.names(colnames(chipseq_12k), unique = TRUE)

#X <- list(rna_seq = rnaseq_all, chipseq = data.matrix(chipseq_12k))
#tokeep <- list(rna_seq = c(10,10), chipseq = c(10,10))


# keep only high variance variables:
# (and sort by variance)
rna_var = sort(apply(rnaseq_all, FUN=var, MARGIN=2), decreasing = T)
rna_var <- rna_var[1:5000]
chipseq_var = sort(apply(chipseq_12k, FUN=var, MARGIN=2), decreasing = T)
chipseq_var <- chipseq_var[1:5000]

rna_seq_best <- rnaseq_all[,names(rna_var)]
chip_seq_best <- chipseq_12k[,names(chipseq_var)]

X <- list(rna_seq = rna_seq_best, chipseq = data.matrix(chip_seq_best))

# one may change number of genes.
tokeep <- list(rna_seq = c(30,30), chipseq = c(30,30))

# run DIABLO:
result.dialbo <- block.splsda(X, y, near.zero.var = TRUE, keepX = tokeep)

# save figures as png
png(file="plots/DIABLO_ADvsCTRL_INDIV.png")

plotIndiv(result.dialbo,
          ind.names = FALSE,
          legend = TRUE, cex = c(1,2),
          title = "AD vs. control")

dev.off()

png(file="plots/DIABLO_ADvsCTRL_VAR.png")
plotVar(result.dialbo, legend = TRUE)
dev.off()

png(file="plots/DIABLO_ADvsCTRL_CIR.png")
circosPlot(result.dialbo, cutoff = 0.85)
dev.off()

corMat <- circosPlot(result.dialbo, cutoff = 0.85)
write.csv(corMat,file="results/corMat_RNAvsCHIP.csv")
