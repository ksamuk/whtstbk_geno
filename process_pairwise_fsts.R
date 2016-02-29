# process a folder of fst log files (from vcftools --weir-pop operations)

# libraries
library("dplyr")
library("stringr")
library("ggplot2")

list.files("functions", full.names = TRUE) %>% sapply(.,source, verbose = FALSE, echo = FALSE) %>% invisible

# list of fst log files

fst_folder <- "data/stats/fst_logs/"
fst_log_files <- list.files(fst_folder)

fst_log_files <- fst_log_files[grep("DK", invert = TRUE, fst_log_files)]

process_fst_log_file <- function(x){
  #print(x)
  file_name <- paste0(fst_folder, x)
  file_split <- strsplit(x, "_") %>% unlist
  
  popA <- file_split %>% .[1:2] %>% paste(collapse = "_")
  popB <- file_split %>% .[3:4] %>% paste(collapse = "_") 
  
  pop1<- sort(c(popA, popB))[1]
  pop2<- sort(c(popA, popB))[2]
  
  type <- file_split %>% .[c(2,4)] %>% sort(decreasing = TRUE) %>% paste(collapse = "_") 
  
  if((file_split[1] == file_split[3])){
    geography <- "sympatric"
  } else{
    geography <- "allopatric"
  }
  
  fst_lines <- scan(file_name, what="", sep="\n", quiet = TRUE)
  
  if (length(fst_lines) == 18){
    fst <- grep("^Weir and Cockerham weighted", fst_lines, value = TRUE) %>% gsub("[^0-9.-]", "", .) %>% as.numeric
    n_inds <- grep("^After filtering.* Individuals", fst_lines, value = TRUE) %>%
      gsub("out of [0-9]* Individuals", "", .) %>%
      gsub("[^0-9.-]", "", .) %>% as.numeric
  } else{
    fst <- NA
    n_inds <- NA
  }
  
  
  data.frame(pop1, pop2, type, geography, fst, n_inds)
  
}

fst_df <- lapply(fst_log_files, process_fst_log_file) %>% rbind_all()

fst_df <- fst_df %>%
  filter(!is.na(fst)) %>%
  filter(!is.nan(fst)) %>%
  filter(pop1 != pop2) %>%
  distinct()

fst_df$fst[fst_df$fst < 0] <- 0

write.table(fst_df, "data/stats/fst_pairwise.txt", row.names = FALSE, quote = FALSE)

