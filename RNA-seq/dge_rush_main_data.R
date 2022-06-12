
library(readr)
library(stringr)

rush_data<- read.csv("/home/rosa/rush_data/all_counts.csv",  header=TRUE)
row.names(rush_data)<- rush_data$hgnc_symbol
rush_data<- rush_data[,2:ncol(rush_data)]

rush_data <- rush_data %>% distinct(hgnc_symbol, .keep_all=TRUE )


rnaseqdata <- read.table("/home/rosa/Dokumente/GSE153873_summary_count.star.txt", sep = "\t",header = TRUE,row.names = 1)
colnames(rnaseqdata) <- paste0(data.frame(str_split(colnames(rnaseqdata),'[.]'))[3,],data.frame(str_split(colnames(rnaseqdata),'[.]'))[1,])
#rnaseq_ercc <- read.table("RNAseq/GSE153873_summary_count.ercc.txt", sep = "\t")
dim(rnaseqdata)
# 27135 different genes and 30 patients


# Groups
AD <- select(rnaseqdata, starts_with("AD"))
young <- select(rnaseqdata, starts_with("Young"))
old <- select(rnaseqdata, starts_with("Old"))


old$hgnc_symbol <- rownames(old)

# 23603 common measured genes 

control_female<- inner_join(rush_data, old, by="hgnc_symbol")