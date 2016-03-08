#install.packages("diveRsity")
library("adegenet")
library("hierfstat")
library("data.table")
library("dplyr")
library("stringr")
library("diveRsity")

### R code
# load diveRsity

data(Test_data, package = "diveRsity")

tmp <- readGenepop(infile = "data/hierfstat/splits_tree_thin.gen")
write.table(Test_data, file = "test.txt", quote = FALSE, col.names = FALSE, row.names = FALSE)

# calculate stats
diff_stats <- diffCalc("data/hierfstat/splits_tree_thin.gen", outfile = "myresults", 
                       fst = TRUE, pairwise = TRUE)
### END