# parse a treemix f4 file

library(dplyr)
library(stringr)

f_file_name <- "data/treemix/fourpop_k10.txt"

#con <- file(description = f_file, open="r")
f_file <- scan(f_file_name , what="", sep="\n")

# harvest block number and size
n_blocks <- str_extract(f_file[2], "[0-9]* blocks") %>% gsub("[^0-9]", "", .) %>% as.numeric
block_size <- str_extract(f_file[2], "size [0-9]*") %>% gsub("[^0-9]", "", .) %>% as.numeric

f_file <- f_file[-grep("^Estimating", f_file)]
f_file <- f_file[-grep("^total", f_file)]
f_file <- f_file[-grep("^npop\\:", f_file)]

