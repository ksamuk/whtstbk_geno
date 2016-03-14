# per locus fst exploration

################################################################################
# initials
################################################################################

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


################################################################################
# FDR correction
################################################################################

fdr_df <- perm_df %>%
  select(matches("pval"))

# little wrapper for fdrtool
apply_fdr_control <- function(x){
  
  # fdr tool doesn't handle NAs for some reason? killing me here.
  pvals <- x[!is.na(x)]
  
  fdr <- fdrtool(pvals, statistic = "pvalue", plot = FALSE, verbose = FALSE)
  
  x[!is.na(x)] <- fdr$qval
  
  x
  
}


fdr_df <- lapply(fdr_df, apply_fdr_control) %>% as.data.frame %>% data.frame(perm_df[,1:2], .)
names(fdr_df) <- names(fdr_df) %>% gsub("pval", "qval", .)

# apply cutoff
out_df <- data.frame(perm_df[,1:2], fdr_df[,3:8] <= 0.05)

################################################################################
# Plot tiled windows on chromosomes
################################################################################

out_long <- gather(out_df, key = comparison, value = outlier, -chr, -pos)


out_long %>%
  filter(chr %in% c(1,4,7,8, 22)) %>%
  filter(grepl("xtx", comparison)) %>%
  filter(outlier == TRUE) %>%
  mutate(comparison_num = as.numeric(as.factor(comparison)) * 2) %>%
  ggplot(aes(x = pos, xmin = pos, xmax = pos + 74999, ymin = comparison_num-1, ymax = comparison_num, fill = comparison))+
  geom_rect()+
  facet_grid(chr~.)
  

