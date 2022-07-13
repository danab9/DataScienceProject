library(readxl)
library(devtools)
library(ggVennDiagram)
library(ggplot2)
sheets <- excel_sheets("RNAseq/GSE153873_AD.vs.Old_diff.genes.xlsx")
resultupreg <- read_excel("RNAseq/GSE153873_AD.vs.Old_diff.genes.xlsx",sheet = sheets[1])
resultdownreg <- read_excel("RNAseq/GSE153873_AD.vs.Old_diff.genes.xlsx",sheet = sheets[2])
upregulated <- read.table("upregulated.txt",header=TRUE)
downregulated <- read.table("downregulated.txt",header = TRUE)

# overlap of up-/downregulated genes in own results and papers results  
length(intersect(upregulated$genes,unlist(resultupreg[,1])))
#314 of 421
length(intersect(downregulated$genes,unlist(resultdownreg[,1])))
#344 of 434


# a lot of overlaps but also some up/downregulated genes but also a lot of genes missing
# in the 50 most significant downregulated genes: overlap of 49 genes
# in the 50 most significant upregulated genes: overlap of 45 genes



x <- list(paper=unlist(resultupreg[,1]),own=upregulated$genes)

venn <- Venn(x)
data <- process_data(venn)
png(file = "VennUpregualted.png", height = 400,width=600)
ggplot() +
  # 1. region count layer
  geom_sf(aes(fill = count), data = venn_region(data),size=2) +
  # 2. set edge layer
  geom_sf(aes(color = name), data = venn_setedge(data), show.legend = TRUE, size = 2) +
  # 3. set label layer
  geom_sf_text(aes(label = name), data = venn_setlabel(data),size=5) +
  # 4. region label layer
  geom_sf_label(aes(label = paste0(count, " (", scales::percent(count/sum(count), accuracy = 2), ")")), 
                data = venn_region(data),
                size = 5) +
   scale_fill_gradient(low = "grey80", high = "red")+
   scale_color_manual(values = c("own" = "red","paper" ="grey70"))+
  theme_void()
dev.off()

x <- list(paper=unlist(resultdownreg[,1]),own=downregulated$genes)

venn <- Venn(x)
data <- process_data(venn)
png(file = "VennDownregulated.png", height = 400,width=600)
ggplot() +
  # 1. region count layer
  geom_sf(aes(fill = count), data = venn_region(data),size=2) +
  # 2. set edge layer
  geom_sf(aes(color = name), data = venn_setedge(data), show.legend = TRUE, size = 2) +
  # 3. set label layer
  geom_sf_text(aes(label = name), data = venn_setlabel(data),size=5) +
  # 4. region label layer
  geom_sf_label(aes(label = paste0(count, " (", scales::percent(count/sum(count), accuracy = 2), ")")), 
                data = venn_region(data),
                size = 5) +
  scale_fill_gradient(low = "grey80", high = "red")+
  scale_color_manual(values = c("own" = "red","paper" ="grey70"))+
  theme_void()
dev.off()

#GO enrichment
paperGO <- read.table('go_genesPaper.txt')
ownGO <- read.table('go_genes.txt')

x <- list(paper=unlist(paperGO),own=unlist(ownGO))

venn <- Venn(x)
data <- process_data(venn)
png(file = "GOenrichment.png", height = 400,width=600)
ggplot() +
  # 1. region count layer
  geom_sf(aes(fill = count), data = venn_region(data),size=2) +
  # 2. set edge layer
  geom_sf(aes(color = name), data = venn_setedge(data), show.legend = TRUE, size = 2) +
  # 3. set label layer
  geom_sf_text(aes(label = name), data = venn_setlabel(data),size=5) +
  # 4. region label layer
  geom_sf_label(aes(label = paste0(count, " (", scales::percent(count/sum(count), accuracy = 2), ")")), 
                data = venn_region(data),
                size = 5) +
  scale_fill_gradient(low = "grey80", high = "red")+
  scale_color_manual(values = c("own" = "red","paper" ="grey70"))+
  theme_void()
dev.off()
