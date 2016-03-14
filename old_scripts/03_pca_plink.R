# plot PCA of genotype data

library(adegenet)
library(pegas)
library(dplyr)
library(stringr)
library(ggplot2)
library(ggthemes)

# a plink 'raw' file (--recodeA
plink_file <- list.files("data/other_formats", pattern = "thinned.raw$", full.names = TRUE)

# read into a genlight object
wht_gen <- read.PLINK(plink_file, parallel = FALSE)

# perform PCA
wht_pca <- glPca(wht_gen, parallel = FALSE, nf = 8)

# extract scores
pca_df <- wht_pca$scores %>% data.frame
rownames(pca_df) <- NULL

# grab individual games
pca_df$id <- rownames(wht_pca$scores)

# prettify names
pca_df$id <- pca_df$id %>% 
  gsub("whtstbk_gbs_2012_brds_", "2012_", .) %>%
  gsub("\\.", "_", .) %>%
  str_trim

# make pop categories
pca_df$pop <- pca_df$id %>%
  gsub("[^A-Z]*", "", .)

# make regions
pop_unique <- pca_df$pop %>% unique
pop_to_region <- c("AN", "CB", "CB", "CB", "AN", "AN", "CB", "CB", "CB", "GY", "GY", "CB",
                   "GY", "GY", "CB", "GY", "GY", "CB", "GY", "GY", "AN", "GY", "AN", "AN")
pca_df$region <- pop_to_region[match(pca_df$pop, pop_unique)]

# add incomplete sex info
sex_df <- read.csv("metadata/whtcmn_sex.csv", stringsAsFactors = FALSE)
sex_df <- sex_df %>%
  mutate(id = paste0(population, individual)) %>%
  select(id, sex)

pca_df <- left_join(pca_df, sex_df)

pca_df$sex[grepl("2012", pca_df$id)] <- "M"

# write metadata file
meta <- pca_df %>%
  select(id, pop, region, sex)

write.csv(meta, file = "metadata/sex_reg.csv", row.names = FALSE)

# rearrange df


# plot PCA
pca_df %>%
  ggplot(aes(x = PC1, y = PC2, color = sex, label = id)) +
  geom_text(size = 3)
  #geom_point(size = 3)

tmp <- dapc.genlight(wht_gen[!is.na(pca_df$sex)], pop = pca_df$sex[!is.na(pca_df$sex)], parallel = FALSE)

