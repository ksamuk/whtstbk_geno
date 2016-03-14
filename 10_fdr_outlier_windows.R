# per locus fst exploration

################################################################################
# initials
################################################################################

library("qqman")
library("dplyr")
library("ggplot2")
library("tidyr")
library("ggthemes")
library("fdrtool")

list.files("functions", full.names = TRUE) %>% sapply(.,source, verbose = FALSE, echo = FALSE) %>% invisible

# the raw "stats" files
stats_folder <- "data/stats/2014"
stats_files <- list.files(stats_folder, pattern = ".txt", full.names = TRUE)

################################################################################
# read in raw data
################################################################################

# master (snp) stats file
fx_df <- read.table("data/stats/snp_stats_master.txt", header = TRUE, stringsAsFactors = FALSE)

# permutation pvals 
perm_df <- read.table("data/stats/outlier_permutation_df.txt", header = TRUE, stringsAsFactors = FALSE)

fdr <- fdrtool(perm_df$fst_wht_cmn_pval, statistic = "pvalue")

fdr


