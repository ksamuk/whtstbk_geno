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

# master (window) stats file
stat_df <- read.table("data/stats/75k_stats_master.txt", header = TRUE, stringsAsFactors = FALSE)
snp_df <- read.table("data/stats/snp_stats_master.txt", header = TRUE, stringsAsFactors = FALSE)


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

out_df$fst_outlier_adj <- with(out_df,(fst_wht_cmn_qval|fst_wht_cbr_qval)&(!fst_cbr_cmn_qval))
out_df$fst_outlier_super <- with(out_df,(fst_wht_cmn_qval&fst_wht_cbr_qval)&(!fst_cbr_cmn_qval))
out_df$xtx_outlier_adj <- with(out_df,(xtx_wht_cmn_qval|xtx_wht_cbr_qval)&(!xtx_cmn_cbr_qval))
out_df$xtx_outlier_super <- with(out_df,(xtx_wht_cmn_qval&xtx_wht_cbr_qval)&(!xtx_cmn_cbr_qval))
out_df$outlier_ultra <- with(out_df,fst_outlier_super & xtx_outlier_super)

out_long <- gather(out_df, key = comparison, value = outlier, -chr, -pos)


out_long %>%
  #filter(chr %in% c(1,4,7,8, 22)) %>%
  filter(grepl("super", comparison)) %>%
  filter(outlier == TRUE) %>%
  mutate(comparison_num = as.numeric(as.factor(comparison)) * 2) %>%
  ggplot(aes(x = pos, xmin = pos, xmax = pos + 74999, ymin = comparison_num-1, ymax = comparison_num, fill = comparison))+
  geom_rect()+
  facet_grid(chr~.)+
  xlim(0, NA)
  
################################################################################
# Join stats and permuted outlier windows
################################################################################# 

# force name compatibility
names(out_df)[2] <- "pos1"
stat_df$chr <- stat_df$chr %>% gsub("chr", "", .) %>% gsub("Un", "XXII", .) %>% as.roman %>% as.numeric

sel_df <- left_join(out_df %>% select(-matches("qval")), stat_df %>% select(-pos2))
sel_long <- gather(sel_df, key = stat, value = value, -chr, -pos1, -matches("adj|super|ultra"))

sel_long %>% 
  ggplot(aes(x = outlier_ultra, y = value, color = stat))+
  geom_jitter()+
  facet_wrap(~stat, scales = "free_y")

# compare ihs/ies/rsb/load/r2 SNP stats for windows

snp_df$pos2 <- (((snp_df$pos / 75000) %>% floor) + 1)*75000
snp_df$pos1 <- snp_df$pos2 - 74999

snp_df <- left_join(out_df %>% select(-matches("qval")), snp_df %>% select(-pos2))
snp_long <- gather(snp_df , key = stat, value = value, -chr, -pos1, -matches("adj|super|ultra"))

snp_long %>% 
  ggplot(aes(x = outlier_ultra, y = value, color = stat))+
  geom_boxplot()+
  facet_wrap(~stat, scales = "free_y")+
  theme(legend.position ="none")


