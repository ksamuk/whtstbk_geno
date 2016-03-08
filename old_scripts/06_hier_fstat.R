#library("adegenet")
library("hierfstat")
#library("data.table")
library("dplyr")
library("stringr")
library("ggplot2")
library("tidyr")

list.files("functions", full.names = TRUE) %>% sapply(.,source, verbose = FALSE, echo = FALSE) %>% invisible

# load hierfstat data
# 'wht_fstat'
load("data/adegenet/wht_bi_snps_hierfstat.R")

# read in metadata
meta_df <- read.csv("metadata/mega_meta.csv")

# save original populations
pop_orig <- wht_fstat$pop

#try: replacing pop with species?

#wht_fstat$pop <- meta_df$species

# try: regions

GY <- c("CL", "MH", "SF", "SR", "SH")
CB <- c("BR", "GC", "LN", "MR", "RT", "SK")

wht_fstat <- wht_fstat %>%
  filter(pop %in% c(GY, CB))

meta_df <- meta_df %>%
  filter(pop %in% c(GY, CB))

meta_df$pop == wht_fstat$pop

wht_fstat$pop <- ifelse(wht_fstat$pop %in% GY, "GY", "CB")
wht_fstat$species <- meta_df$species
wht_fstat$year <- meta_df$year
wht_fstat$sex <- meta_df$gen.sex

wht_fstat <- wht_fstat %>%
  filter(year != 2012) %>%
  select(-year)

wht_fstat_GY <- wht_fstat %>%
  filter(pop == "GY") %>%
  mutate(pop = species) %>%
  select(-species, -sex)

wht_fstat_GY_male <- wht_fstat %>%
  filter(pop == "GY", sex == "M") %>%
  mutate(pop = species) %>%
  select(-species, -sex)

wht_fstat_GY_female <- wht_fstat %>%
  filter(pop == "GY", sex == "F") %>%
  mutate(pop = species) %>%
  select(-species, -sex)
  
wht_fstat_CB <- wht_fstat %>%
  filter(pop == "CB") %>%
  mutate(pop = species) %>%
  select(-species, -sex)

wht_fstat_CB_male <- wht_fstat %>%
  filter(pop == "CB", sex == "M") %>%
  mutate(pop = species) %>%
  select(-species, -sex)

wht_fstat_CB_female <- wht_fstat %>%
  filter(pop == "CB", sex == "F") %>%
  mutate(pop = species) %>%
  select(-species, -sex)

# calculate fst per locus
fst_GY_male <- fst_per_locus(wht_fstat_GY_male, "cmn", "wht")
fst_GY_female <- fst_per_locus(wht_fstat_GY_female, "cmn", "wht")
fst_GY_both <- fst_per_locus(wht_fstat_GY, "cmn", "wht")

fst_CB_male <- fst_per_locus(wht_fstat_CB_male, "cmn", "wht")
fst_CB_female <- fst_per_locus(wht_fstat_CB_female, "cmn", "wht")
fst_CB_both <- fst_per_locus(wht_fstat_CB, "cmn", "wht")

write.csv(fst_GY_male, "data/hierfstat/fst_GY_male.csv", row.names = FALSE, quote = FALSE)
write.csv(fst_GY_female, "data/hierfstat/fst_GY_female.csv", row.names = FALSE, quote = FALSE)
write.csv(fst_GY_both, "data/hierfstat/fst_GY_both.csv", row.names = FALSE, quote = FALSE)

write.csv(fst_CB_male, "data/hierfstat/fst_CB_male.csv", row.names = FALSE, quote = FALSE)
write.csv(fst_CB_female, "data/hierfstat/fst_CB_female.csv", row.names = FALSE, quote = FALSE)
write.csv(fst_CB_both, "data/hierfstat/fst_CB_both.csv", row.names = FALSE, quote = FALSE)

fst_GY_male <- read.csv("data/hierfstat/fst_GY_male.csv")
fst_GY_female <-read.csv("data/hierfstat/fst_GY_female.csv")
fst_GY_both <-read.csv("data/hierfstat/fst_GY_both.csv")

fst_CB_male <-read.csv("data/hierfstat/fst_CB_male.csv")
fst_CB_female <-read.csv("data/hierfstat/fst_CB_female.csv")
fst_CB_both <-read.csv("data/hierfstat/fst_CB_both.csv")

fst_GY_male$fst[fst_GY_male$fst < 0] <- 0
fst_GY_female$fst[fst_GY_female$fst < 0] <- 0
fst_GY_both$fst[fst_GY_both$fst < 0] <- 0

fst_CB_male$fst[fst_CB_male$fst < 0] <- 0
fst_CB_female$fst[fst_CB_female$fst < 0] <- 0
fst_CB_both$fst[fst_CB_both$fst < 0] <- 0


# plot results

fst_GY_male$fst.outlier <- is.outlier(fst_GY_male$fst, cutoff = 0.99)
fst_GY_female$fst.outlier<- is.outlier(fst_GY_female$fst, cutoff = 0.99)
fst_GY_both$fst.outlier<- is.outlier(fst_GY_both$fst, cutoff = 0.99)

fst_CB_male$fst.outlier<- is.outlier(fst_CB_male$fst, cutoff = 0.99)
fst_CB_female$fst.outlier<- is.outlier(fst_CB_female$fst, cutoff = 0.99)
fst_CB_both$fst.outlier<- is.outlier(fst_CB_both$fst, cutoff = 0.99)


fst_combined <- data.frame(chr = fst_CB_both$chr, pos = fst_CB_both$pos, 
                           fst_GY_male = fst_GY_male$fst, fst_GY_male_outlier = fst_GY_male$fst.outlier,
                           fst_GY_female = fst_GY_female$fst, fst_GY_female_outlier = fst_GY_female$fst.outlier,
                           fst_CB_male = fst_CB_male$fst, fst_CB_male_outlier = fst_CB_male$fst.outlier,
                           fst_CB_female = fst_CB_female$fst, fst_CB_female_outlier = fst_CB_female$fst.outlier)

fst_combined_long1 <- gather(fst_combined, key = fst_class, value = fst, fst_GY_male, fst_GY_female, fst_CB_male, fst_CB_female)
fst_combined_long2 <- gather(fst_combined, key = fst_outlier_class, value = fst_outlier, fst_GY_male_outlier, fst_GY_female_outlier, fst_CB_male_outlier, fst_CB_female_outlier)

fst_combined_long <- data.frame(chr = fst_combined_long1$chr, pos = fst_combined_long1$pos,
                           fst_class = fst_combined_long1$fst_class, fst = fst_combined_long1$fst,
                           fst_outlier_class = fst_combined_long2$fst_outlier_class, fst_outlier = fst_combined_long2$fst_outlier)


fst_combined_long %>%
  sample_frac(1) %>%
  #filter(chr == 1) %>%
  #filter(fst.outlier.male == TRUE | fst.outlier.female == TRUE) %>%
  ggplot(aes(x = pos, y = fst, color = fst_outlier)) +
  geom_point()+
  #stat_smooth(span = 0.2, n = 100,  se = FALSE)+
  facet_grid(fst_class~chr)

fst_combined_long %>%
  sample_frac(1) %>%
  filter(chr == 11) %>%
  filter(fst_outlier == TRUE, fst_class == "fst_GY_male") %>%
  arrange(pos)
  #filter(fst.outlier.male == TRUE | fst.outlier.female == TRUE) %>%
  ggplot(aes(x = pos, y = as.numeric(fst_outlier), color = fst_outlier_class)) +
  stat_smooth(n = 1000, span = 0.1, se = FALSE)+
  facet_wrap(~chr)

fst_combined %>% 
  filter(fst_CB_male > 0.3)

head(fst_whtcmn)
