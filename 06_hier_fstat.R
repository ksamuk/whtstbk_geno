#library("adegenet")
library("hierfstat")
library("data.table")
library("dplyr")
library("stringr")
library("ggplot2")

list.files("functions", full.names = TRUE) %>% sapply(.,source, verbose = FALSE, echo = FALSE) %>% invisible

# load hierfstat data
# 'wht_fstat'
load("data/adegenet/wht_bi_snps_hierfstat.R")

# read in metadata
meta_df <- read.csv("metadata/mega_meta.csv")

# save original populations
pop_orig <- wht_fstat$pop

#try: replacing pop with species?
wht_fstat$pop <- meta_df$species

#do by...sex?
wht_fstat$sex <- meta_df$gen.sex

wht_fstat_male <- wht_fstat %>%
  filter(sex == "M") %>%
  select(-sex)

wht_fstat_female <- wht_fstat %>%
  filter(sex == "F") %>%
  select(-sex)

wht_fstat <- wht_fstat %>%
  select(-sex)

# calculate fst per locus
fst_whtcmn_male <- fst_per_locus(wht_fstat_male, "cmn", "wht")
fst_whtcmn_female <- fst_per_locus(wht_fstat_female, "cmn", "wht")
fst_whtcmn <- fst_per_locus(wht_fstat, "cmn", "wht")

write.csv(fst_whtcmn_male, "data/hierfstat/wht_cmn_fst_male.csv", row.names = FALSE, quote = FALSE)
write.csv(fst_whtcmn_female, "data/hierfstat/wht_cmn_fst_female.csv", row.names = FALSE, quote = FALSE)
write.csv(fst_whtcmn, "data/hierfstat/wht_cmn_fst.csv", row.names = FALSE, quote = FALSE)

fst_whtcmn_male  <- read.csv("data/hierfstat/wht_cmn_fst_male.csv", stringsAsFactors = FALSE)
fst_whtcmn_female <- read.csv("data/hierfstat/wht_cmn_fst_female.csv", stringsAsFactors = FALSE)
fst_whtcmn <- read.csv("data/hierfstat/wht_cmn_fst.csv", stringsAsFactors = FALSE)

# plot results

fst_whtcmn_male$fst.outlier <- is.outlier(fst_whtcmn_male$fst)
fst_whtcmn_female$fst.outlier <- is.outlier(fst_whtcmn_female$fst)

fst_combined1 <- fst_whtcmn_male %>%
  select(chr, pos, fst, fst.outlier) %>%
  rename(fst.outlier.male = fst.outlier, fst.male = fst)

fst_combined2 <- fst_whtcmn_female %>%
  select(chr, pos, fst, fst.outlier) %>%
  rename(fst.outlier.female = fst.outlier, fst.female = fst)

fst_combined <- left_join(fst_combined1, fst_combined2)


fst_combined %>%
  sample_frac(0.1) %>%
  ggplot(aes(x = pos, y = fst.male)) +
  geom_point()+
  geom_point(aes(x = pos, y = fst.female))+
  facet_wrap(~chr)

fst_whtcmn_female %>%
  ggplot(aes(x = fst)) +
  geom_histogram()+
  facet_wrap(~chr)




