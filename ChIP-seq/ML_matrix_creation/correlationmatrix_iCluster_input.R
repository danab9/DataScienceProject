#get correlation matrix of all ranges
file <- read.table("C:/Users/georg/Documents/Uni/Semester_2/Data Science in Life Sciences/Main_project/iCluster_prep/second_Laptop/iCluster_mat_input.txt")
install.packages("psych")
library(psych)

mat_cor <- data.matrix(file[,4:30])
corPlot(mat_cor)

corPlot(file[,4:30])
