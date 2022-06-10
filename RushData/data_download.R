library(readr)

report<- read.csv("/home/rosa/Downloads/experiment_report_2022_6_10_10h_48m.tsv", sep="\t", header=TRUE)
samples<- read.table("/home/rosa/rush_data/sample.txt", sep="\t", header=TRUE)
filenames<- read.table("/home/rosa/rush_data/filenames.txt")


report["file_name"]<- NA
for (filename in filenames$V1){
  for (i in 1:nrow(report)){
    list<- report[i,15]
    if (grepl(filename,list )){
      report[i, ncol(report)] <- filename
      
    }
  }
  
  
}