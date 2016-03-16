# per locus fst exploration

################################################################################
# initials
################################################################################

library("dplyr")
library("ggplot2")
library("tidyr")
library("ggthemes")
library("fdrtool")
library("sinaplot")

list.files("functions", full.names = TRUE) %>% sapply(.,source, verbose = FALSE, echo = FALSE) %>% invisible

# the raw "stats" files
stats_folder <- "data/stats/2014"
stats_files <- list.files(stats_folder, pattern = ".txt", full.names = TRUE)

################################################################################
# read in raw data
################################################################################

# master (window) stats file
wind_df <- read.table("data/stats/75k_stats_master.txt", header = TRUE, stringsAsFactors = FALSE)
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
out_df$fst_outlier_triple <- with(out_df,(fst_wht_cmn_qval|fst_wht_cbr_qval)&(!fst_cbr_cmn_qval))
out_df$xtx_outlier_adj <- with(out_df,(xtx_wht_cmn_qval|xtx_wht_cbr_qval)&(!xtx_cmn_cbr_qval))

out_long <- gather(out_df, key = comparison, value = outlier, -chr, -pos)

out_long %>%
  #filter(chr %in% c(1,4,7,8, 22)) %>%
  filter(grepl("adj", comparison)) %>%
  filter(outlier == TRUE) %>%
  mutate(comparison_num = as.numeric(as.factor(comparison)) * 2) %>%
  ggplot(aes(x = pos, xmin = pos, xmax = pos + 74999, ymin = comparison_num, ymax = comparison_num-2, fill = comparison))+
  geom_rect()+
  facet_grid(chr~.)+
  xlim(0, NA)+
  theme(legend.position ="none")
  
################################################################################
# Join stats and permuted outlier windows
################################################################################# 

# force name compatibility
names(out_df)[2] <- "pos1"
wind_df$chr <- wind_df$chr %>% gsub("chr", "", .) %>% gsub("Un", "XXII", .) %>% as.roman %>% as.numeric

sel_df <- left_join(out_df %>% select(chr, pos1, matches("adj")), wind_df %>% select(-pos2))
sel_long <- gather(sel_df, key = stat, value = value, -chr, -pos1, -matches("adj"))

sel_long %>% 
  ggplot(aes(x = fst_outlier_adj, y = value, fill = stat))+
  geom_violin()+
  facet_wrap(~stat, scales = "free_y")

# compare ihs/ies/rsb/load/r2 SNP stats for windows

snp_df <- read.table("data/stats/snp_stats_master.txt", header = TRUE, stringsAsFactors = FALSE)

snp_df$pos2 <- (((snp_df$pos / 75000) %>% floor) + 1)*75000
snp_df$pos1 <- snp_df$pos2 - 74999

#snp_df <- left_join(out_df %>% select(-matches("qval")), snp_df %>% select(-pos2, -matches("xtx|fst|outlier|pval")))
snp_df <- left_join(out_df %>% select(-matches("qval")), snp_df %>% select(-pos2))


# pi site data
wht_pi <- read.table("data/stats/wht_pi.site.2014.txt", header = TRUE, stringsAsFactors = FALSE)
cmn_pi <- read.table("data/stats/cmn_pi.site.2014.txt", header = TRUE, stringsAsFactors = FALSE)
cbr_pi <- read.table("data/stats/cbr_pi.site.2014.txt", header = TRUE, stringsAsFactors = FALSE)

pi_df <- left_join(wht_pi, cmn_pi, by = c("CHROM", "POS"))
pi_df <- left_join(pi_df, cbr_pi, by = c("CHROM", "POS"))

names(pi_df) <- c("chr", "pos", "pi_site_wht", "pi_site_cmn", "pi_site_cbr")

pi_df$chr <- pi_df$chr %>% gsub("chr", "", .) %>% gsub("Un", "XXII", .) %>% as.roman %>% as.numeric

snp_df <- left_join(snp_df, pi_df)

################################################################################
# Synthetic variables
################################################################################# 

snp_df$load_abs_ev1 <- abs(snp_df$load_ev1)
snp_df$load_abs_ev2 <- abs(snp_df$load_ev2)
snp_df$r2_delta <- snp_df$r2_wht - snp_df$r2_cmn


# ihs delta/abs stats (not that informative, largely superceded by rsb, i think)
snp_df$ihs_delta <- snp_df$ihs_wht - snp_df$ihs_cmn
snp_df$ihs_delta_abs <- abs(snp_df$ihs_wht - snp_df$ihs_cmn)
snp_df$ihs_abs_wht <- abs(snp_df$ihs_wht)
snp_df$ihs_abs_cmn <- abs(snp_df$ihs_cmn)
snp_df$ihs_abs_delta <- snp_df$ihs_abs_wht - snp_df$ihs_abs_cmn


# join in windowed stats
wind_df$chr <- wind_df$chr <- wind_df$chr %>% gsub("chr", "", .) %>% gsub("Un", "XXII", .) %>% as.roman %>% as.numeric

snp_df <- left_join(snp_df, wind_df %>% select(-matches("fst|r2")))

snp_long <- gather(snp_df , key = stat, value = value, -chr, -pos1, -pos, -matches("adj"))

################################################################################
# plots
################################################################################# 

# whole chromosome view
snp_long %>% 
  filter(chr == 4) %>%
  filter(!grepl("cbr", stat)) %>%
  filter(grepl("rsb_wht_cmn|fst_wht_cmn|fst_cbr_cmn|tajd_wht|tajd_cmn|pi_site_wht|pi_site_cmn", stat)) %>%
  ggplot(aes(x = pos, y = value, color = xtx_outlier_adj))+
  geom_point(alpha = 1)+
  facet_grid(stat~., scales = "free_y")+
  theme_classic()+
  theme(axis.text.x = element_blank(), 
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        panel.margin = unit(0.1, "cm"),
        strip.background = element_blank(),
        legend.position ="none")+
  scale_color_manual(values=c("grey", "red"))

# boxplot outlier vs. genome
snp_long %>% 
  #filter(chr == 21) %>%
  filter(!grepl("cbr", stat)) %>%
  filter(grepl("rsb_wht_cmn|fst_wht_cmn|fst_cbr_cmn|tajd_wht|tajd_cmn|pi_site_wht|pi_site_cmn", stat)) %>%
  mutate(outlier = fst_outlier_adj|xtx_outlier_adj) %>%
  ggplot(aes(x = outlier, y = value, color = stat))+
  geom_boxplot()+
  #geom_jitter(alpha = 0.5)+
  facet_wrap(~stat, scales = "free_y")+
  theme_classic()+
  theme(#axis.text.x = element_blank(), 
        axis.ticks.x = element_blank(),
        #axis.title.x = element_blank(),
        panel.margin = unit(0.1, "cm"),
        strip.background = element_blank(),
        legend.position ="none")


snp_long %>%
  #filter(chr %in% c(1,4,7,8, 22)) %>%
  filter(grepl("adj", stat)) %>% 
  filter(outlier == TRUE) %>%
  mutate(comparison_num = as.numeric(as.factor(comparison)) * 2) %>%
  ggplot(aes(x = pos, xmin = pos, xmax = pos + 74999, ymin = comparison_num-1, ymax = comparison_num, fill = comparison))+
  geom_rect()+
  facet_grid(chr~.)+
  xlim(0, NA)


