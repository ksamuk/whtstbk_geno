# inject chromosome names into a plink map file
# vcf tools fails to properly convery chrosome names, hence this script


args <- commandArgs(trailingOnly = TRUE)

map_file_name <- args[1]
#map_file_name <- "temp_plink_pre_prune.map"

map_file <- read.table(map_file_name)
tmp <- gsub(":", ",", map_file[,2])
map_file[,1] <- unlist(lapply(strsplit(tmp, split = ","), function(x)x[[1]]))
map_file[,1] <- gsub("chr", "", map_file[,1])
map_file[,1][map_file[,1] == "Un"] <- "XXII"
map_file[,1] <- as.numeric(as.roman(map_file[,1]))

write.table(map_file, map_file_name, col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")
