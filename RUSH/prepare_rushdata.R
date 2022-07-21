library(readr)
library(dplyr)
library('biomaRt')
library(tximport)

# info about samples rna seq data from rush , filtering out file name of the samples
# the resulting overview tables are in the folder overviews
# the code is locally only runnable. The code  just extract information from the rather big table, that you get when you download files from ENCODE
report<- read.csv("/home/rosa/Downloads/experiment_report_2022_6_10_10h_48m.tsv", sep="\t", header=TRUE)
samples<- read.table("/home/rosa/rush_data/sample.txt", sep="\t", header=TRUE)
filenames<- read.table("~/DataScienceProject/RushData/data/filenames_ADF.txt")


report["file_name"]<- NA
for (filename in filenames$V1){
  for (i in 1:nrow(report)){
    list<- report[i,15]
    if (grepl(filename,list )){
      report[i, ncol(report)] <- filename
      
    }
  }
  
}
reduced_report<- report[, c(1,2,17,40, 3,4,9,21,22,23, 18,19)]

write_csv(reduced_report, "~/DataScienceProject/RushData/overviews/ADfemale_overview.csv")


# in the following , we extract the expected counts from the files and convert ensembl ids to hgnc symbols
# the rsem tsv files are zipped
# dir saves path to the locally stored rsem files

dir<- "/home/rosa/rush_data/ADfemale"
overview<- read.csv("~/DataScienceProject/RushData/overviews/ADfemale_overview.csv", header=TRUE)

# combining all 29 rna seq data to one rsem object

files <- file.path(dir, paste0(overview$file_name, ".tsv"))
names(files) <-overview$Accession
txi.rsem <- tximport(files, type = "rsem", txIn = FALSE, txOut = FALSE)

# expected count table 
rush_counts<- as.data.frame(txi.rsem$counts)

# removing zero counts and ercc counts
rush_counts<- rush_counts[650:nrow(rush_counts),]
rush_counts<- rush_counts[1:58780,]
#removing par_Y counts 
rush_counts<- filter(rush_counts, !endsWith(rownames(rush_counts), "PAR_Y"))

# remove endings of gene names 
row.names(rush_counts)<- gsub("\\..*", "", rownames(rush_counts))

# convert ensembl ids to gene symbols
mart <- useMart(dataset="hsapiens_gene_ensembl","ensembl")
G_list <- getBM(filters= "ensembl_gene_id", 
                attributes=c( 'ensembl_gene_id',
                              'hgnc_symbol', 'gene_biotype'),
                values= rownames(rush_counts),
                mart= mart)

rush_counts$ensembl_gene_id <- rownames(rush_counts)

merged<- left_join(rush_counts, G_list, by="ensembl_gene_id")



# filter out  data set

# all not protein coding without a gene symbol 
not_protein_coding<- merged$ensembl_gene_id[which(merged$gene_biotype!= "protein_coding" & merged$hgnc_symbol == "")]
only_ensembl_info<- merged$ensembl_gene_id[which(is.na(merged$gene_biotype) & is.na(merged$hgnc_symbol))]
merged<- filter(merged, !(ensembl_gene_id %in% c(not_protein_coding, only_ensembl_info)))


no_symbol_yet <- filter(merged, hgnc_symbol =="")

G_list <- getBM(filters= "ensembl_gene_id", 
                attributes=c( 'refseq_mrna', 'ensembl_gene_id',
                              'hgnc_symbol', 'gene_biotype'),
                values= no_symbol_yet$ensembl_gene_id,
                mart= mart)



# remove all rows without gene symbol and no  refseq id
no_info<- G_list$ensembl_gene_id[which(G_list$refseq_mrna == "")]
merged<- filter(merged, !(ensembl_gene_id %in% no_info))



no_symbol_yet<- filter(G_list, refseq_mrna != "")

# obtained gene symbol with help of resfseq id using website: .... 
refseq_genesym<- read.csv("~/DataScienceProject/RushData/data/refseq_genesymbol.tsv", header=FALSE)
colnames(refseq_genesym)<- c("refseq_mrna", "hgnc_symbol")

now_symbol_yet <- left_join(no_symbol_yet, refseq_genesym, by="refseq_mrna")

# remove all rows without symbols

no_symbol <- filter(now_symbol_yet, hgnc_symbol.y == "")$ensembl_gene_id
merged<- nrow(filter(merged, !(ensembl_gene_id %in% no_symbol)))



now_symbol_yet <- filter(now_symbol_yet, hgnc_symbol.y != "")
now_symbol_yet<- dplyr::select(now_symbol_yet, -c( hgnc_symbol.y))
now_symbol_yet<- unique(now_symbol_yet)

# add now gene symbol to the ensembl ids 

positions<- which(merged$ensembl_gene_id %in% now_symbol_yet$ensembl_gene_id)

helper<-data.frame(ensembl_gene_id= merged$ensembl_gene_id[positions]  ,
                   row_num=positions)

now_symbol_yet<- filter(left_join(now_symbol_yet, helper, by="ensembl_gene_id"), !is.na(row_num))
now_symbol_yet<- now_symbol_yet[order(now_symbol_yet$row_num),]

merged[which(merged$ensembl_gene_id %in% now_symbol_yet$ensembl_gene_id), "hgnc_symbol"] <- now_symbol_yet$hgnc_symbol

# now every ensembl id has a gene symbol 

# get rid of ensembl ids with same gene symbol 

occur<- as.data.frame(table(merged$hgnc_symbol,useNA = 'always'))

duplicates_hgnc <- filter( occur, Freq >1)


merged<- filter(merged, !( hgnc_symbol %in% duplicates_hgnc$Var1))



merged<- merged  %>% dplyr::select(ensembl_gene_id, hgnc_symbol, gene_biotype,   everything())


merged[, 4:ncol(merged)]<- round(merged[, 4:ncol(merged)])

write.csv(file="~/DataScienceProject/RushData/data/all_counts_female.csv", merged)







