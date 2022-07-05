rnaseqdata <- read.table("~/Dokumente/GSE153873_summary_count.star.txt", sep = "\t",header = TRUE,row.names = 1)

up.genes <- read.table("~/DataScienceProject/RNA-seq/UpDownreg_genes/upregulatedgenes.txt", sep = "\t")
down.genes <- read.table("~/DataScienceProject/RNA-seq/UpDownreg_genes/downregulationgenes.txt", sep = "\t")


all.genes<- rownames(rnaseqdata)
#background.total<- all.genes[ which(!(all.genes %in% c(up.genes$V1, down.genes$V1)))]


write.csv(all.genes, "background.gens.txt", row.names = FALSE, quote = FALSE)



mouse_gene <- read.table("~/enrichments/mousegenes.txt", header=FALSE)
mouse_gene <- c(mouse_gene$V1)



library(biomaRt)
human <- useMart("ensembl", dataset = "hsapiens_gene_ensembl", host = "https://dec2021.archive.ensembl.org/")
mouse <- useMart("ensembl", dataset = "mmusculus_gene_ensembl", host = "https://dec2021.archive.ensembl.org/")

genesV2 = getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", values = mouse_gene , mart = mouse, attributesL = c("hgnc_symbol"), martL = human, uniqueRows=T)
humanx <- unique(genesV2[, 2])

write.csv(as.data.frame(humanx), "~/enrichments/zinc.gens.txt", row.names = FALSE, quote = FALSE)


zinc.related.genes <- read.csv("~/enrichments/zinc_related_genes.tsv", sep="\t", header = TRUE)

zinc.related.genes<- zinc.related.genes %>%
  mutate(gene_type = case_when(logFC > 0 & adj.P.Val < 0.05 ~ "upregulated",
                               logFC < 0 & adj.P.Val < 0.05 ~ "downregulated",
                               TRUE ~ "not significant"))   
print(table(zinc.related.genes$gene_type))


up.zinc.genes <- filter(zinc.related.genes, gene_type  == "upregulated")$Gene.symbol
down.zinc.genes <- filter(zinc.related.genes, gene_type  == "downregulated")$Gene.symbol



intersect(up.genes$V1, up.zinc.genes)
intersect(down.genes$V1, down.zinc.genes)

intersect(all.genes, humanx )




