#install.packages("dplyr")
library(dplyr)

setwd("~/Desktop/coverage")

#create a list of the files from your target directory
file_list <- list.files(path="~/Desktop/coverage")

data_set <- read.table(file_list[1], header=TRUE, sep="\t")
data_set <- data_set[,c("chrom","mean"), drop=FALSE]
names(data_set)[names(data_set) == "mean"] <- toString(file_list[1])

for (file in file_list[-1]){
  temp <- read.table(file, header=TRUE, sep="\t")
  temp <- temp[,c("chrom","mean"), drop=FALSE]
  names(temp)[names(temp) == "mean"] <- toString(file)
  data_set <- full_join(data_set, temp, by = "chrom")
}

write.csv(data_set,file="coverage.csv")


