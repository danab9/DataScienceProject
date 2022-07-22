library(Biobase)
library(limma)

normalizedRNA <- read.table("normalizedRNA.txt",header=TRUE,sep = ",",row.names = 1)
sheets <- excel_sheets("41588_2020_696_MOESM3_ESM.xlsx")
clinical <- data.frame(read_excel("41588_2020_696_MOESM3_ESM.xlsx",sheet = sheets[1],col_names = TRUE))
colnames(clinical) <- clinical[2,]
clinical <- clinical[-c(1,2),]
rownames(clinical) <- paste0(clinical$`Study Group`,'X',clinical$`Sample ID`)
colnames(clinical)[2] <- 'group'

#Only old and AD
normalizedRNA <- dplyr::select(normalizedRNA, starts_with(c('old',"AD")))

# Exoressionset for linear model
readES = ExpressionSet(data.matrix(normalizedRNA))
pData(readES) = clinical[sampleNames(readES),]
readES

#Relating groups (AD and old) to gene expression variation
mm = model.matrix(~group, data=pData(readES))
f1 = lmFit(readES, mm)
ef1 = eBayes(f1)

t <- topTable(ef1)
sum(p.adjust(ef1$p.value[,2],method = 'BH')<0.05)
#640 significant genes

# top 3 genes
k1 <- split(exprs(readES)[rownames(t)[1],], readES$group)
k2 <- split(exprs(readES)[rownames(t)[2],], readES$group)
k3 <- split(exprs(readES)[rownames(t)[3],], readES$group)

d <- data.frame(x = c(unlist(k1),unlist(k2),unlist(k3)), 
                Class = rep(rep(letters[1:length(k1)],times = sapply(k1,length)),3),
                Gene=rep(rownames(t)[1:3],c(length(readES$group),length(readES$group),length(readES$group))))
d$Class[d$Class=='a'] <- 'old'
d$Class[d$Class=='b'] <- 'AD'

png(file = "lmFit.png", height = 400,width=600)
p <- ggplot(d, aes(x=Class, y=x,fill=Gene)) + 
  geom_boxplot(position=position_dodge(0.8))+
  geom_dotplot(binaxis='y', stackdir='center', 
               position=position_dodge(0.8))+
  #geom_dotplot(binaxis='y', stackdir='center', dotsize=1)+
  theme_grey(base_size = 22)
p+scale_fill_brewer(palette="Greys")+
  labs(title="", y = "Gene expression")
dev.off()
# description: boxplot of the gene expression of the top 3 genes (most significant) 
# in AD and old
# x-axis: AD and old group, y-axis: gene expression
