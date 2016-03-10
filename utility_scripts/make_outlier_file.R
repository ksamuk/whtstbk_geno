# create a basic list of outliers
# xtx outliers in this case

library(dplyr)

# the master stats file
stats_df <- read.table("data/stats/snp_stats_master.txt", header = TRUE, stringsAsFactors = FALSE)

non_outlier_df <- stats_df %>%
  filter(!(outlier_xtx_cmn_cbr|outlier_xtx_wht_cbr|outlier_xtx_wht_cmn)) %>%
  select(chr, pos)

write.table(non_outlier_df, "metadata/non_outlier_sites.txt", row.names = FALSE)

outlier_df <- stats_df %>%
  mutate(outlier = (!(outlier_xtx_cmn_cbr) & (outlier_xtx_wht_cbr|outlier_xtx_wht_cmn))) %>%
  select(chr, pos, outlier)

outlier_df %>%
  filter(chr == 18) %>% View

write.table(outlier_df, "metadata/outlier_test.txt", row.names = FALSE)
