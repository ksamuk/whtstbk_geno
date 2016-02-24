# generate a 'popfile' from a metadata file

library("dplyr")

meta_df <- read.csv("metadata/mega_meta.csv", header = TRUE, stringsAsFactors = FALSE)
head(meta_df)

# identify populations with very low sampling
pop_counts <- meta_df %>% group_by(pop) %>% tally %>% filter(n >= 6)
meta_df <- meta_df %>% 
  filter(pop %in% pop_counts$pop)

# list of suspect samples
# suspect = v. high misssing data and/or dissimilarity outliers with respect to samples from same location
suspect_samples <- c("SK36", "LN30", "SR20", "SR15", "CL64_2", "whtstbk_gbs_2012_brds_WR15", "whtstbk_gbs_2012_brds_SF16", "whtstbk_gbs_2012_brds_PQ11", "whtstbk_gbs_2012_brds_CP22", "whtstbk_gbs_2012_brds_CP8")

# apply filters for popfile
meta_df <- meta_df %>%
  select(id, pop, cluster) %>%
  filter(!(id %in% suspect_samples)) %>%
  filter(!grepl("SRX", id))

# write to file
write.table(meta_df, "metadata/popfile_filtered.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)
