# parse a treemix f4 file

# From Treemix Manual (1.1):
# ...The output is four columns. These are the populations used to calculate the f4 statistic, the f4 statistic, 
# the standard error in the f4 statistic,and the Z-score. 
# For example, the following line shows the f4 statistic where Han is population A,
# Colombian is population B, Sardinian is population C, and Dai is population D:
#  Han,Colombian;Sardinian,Dai -0.00773721 0.000730214 -10.5958


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

split_f_file_line <- function(x){
  
  row_df <- data.frame(t(data.frame(strsplit(x, split = ";|,| "))))
  names(row_df) <- c("pop1", "pop2", "pop3", "pop4", "f4", "f4_se", "f4_zscore") 
  row.names(row_df) <- NULL
  row_df
}

f4_df <- do.call("rbind", lapply(f_file, split_f_file_line))


f4_df_cluster <- f4_df[,1:4] %>% gsub("[A-Z]*|_", "", .)
