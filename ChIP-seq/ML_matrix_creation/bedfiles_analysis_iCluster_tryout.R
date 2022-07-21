#For examination of the merrged bed files

#read in bed file containing all bed fies merged, which will be used for iCLuster method
file <- read.table("/home/georgvonarnim/Dokumente/DataScience/iCluster/multibed.bed")
summary(file)

#create histogram of the score of the bed files
hist(file$V4, main="score distribution of the bed files", xlab="score")

#maximal score
max(file$V4)

#calculate length of the reads
file$V5 <- file$V3 - file$V2

#only plot for up to 2000
hist(file$V5[file$V5<2000])

#longest read
max(file$V5)






### For creation of Matrix, which is required for iClusterBayes

mat <- read.table("/home/georgvonarnim/Dokumente/DataScience/iCluster/iCluster_matrix.bed")

hist(rowSums(mat != 0))
sum(rowSums(mat != 0)>9)




#all in all 27 patients, three rows for position and chromosome -> check for how many samples coverage information is present
sum(rowSums(mat != 0)>17)   

sum(rowSums(mat != 0)>=27)
mat_27 <- mat[which(rowSums(mat != 0)>=17),]
range_length2 = mat_27$V3 - mat_27$V2
hist(range_length2[range_length2 < 15])
sum(range_length2>15)

mat_new <- mat_27[which(range_length2>15),]
write.table(mat_new, "/home/georgvonarnim/Dokumente/DataScience/iCluster/iCluster_matrix_firsttry_100k.bed", col.names=FALSE, row.names=FALSE, sep = "\t")



#get overview over the maximum values 
list <- apply(mat[ , 4:30], 1, max)
summary(rowSums(mat[ , c(4:30)], na.rm=TRUE))
sums_max <- rowSums(mat[ , c(4:30)], na.rm=TRUE)
hist(sums_max[sums_max<100])


#number of rows with a maximum coverage depth of 20 
sum(sums_max<20)
sum(sums_max>20)


#number of rows with a maximum coverage depth of 50 
sum(sums_max<50)

sum(rowSums(mat != 0)>7)   

sum(rowSums(mat != 0)>7 & sums_max>20)




mat_max_10 <- mat[which(rowSums(mat != 0)>10),]

write.table(mat_max_10, "/home/georgvonarnim/Dokumente/DataScience/iCluster/iCluster_matrix_10rows.bed", col.names=FALSE, row.names=FALSE, sep = "\t")



# get overview on how long the ranges are in which the coverage is present (based on the margingform of unionbedg from bedtools)

range_length = mat_max_10$V3 - mat_max_10$V2
hist(range_length[range_length < 15])
sum(range_length<15)

mat_length_15 <- mat_max_10[which(range_length>14),]
#sum(rowSums(mat[ , c(4:30)], na.rm=TRUE) < 10)
#sum(rowSums(mat_max_10 != 0)<10)

write.table(mat_length_15, "/home/georgvonarnim/Dokumente/DataScience/iCluster/iCluster_mat_length_15.bed", col.names=FALSE, row.names=FALSE, sep = "\t")

list <- apply(mat_length_15[ , 4:30], 1, max)  
sum(list<15)



#get correlation matrix of all ranges

install.packages("psych")
library(psych)

mat_cor <- data.matrix(mat[,4:30])
corPlot(mat_cor)

corPlot(data.matrix(mat_27[,4:30]))
