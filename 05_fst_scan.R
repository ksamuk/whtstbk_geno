library("adegenet")
library("hierfstat")
library("data.table")
library("dplyr")
library("stringr")

plink_file <- list.files("data/other_formats", pattern = "full", full.names = TRUE)

# read into a genlight object
wht_gen <- read.PLINK(plink_file, parallel = FALSE)

#meta data df
meta_df <- read.csv("metadata/sex_reg.csv")

pop_file <- meta_df %>%
  select(id, pop)

write.table(pop_file, "metadata/wht_cmn_pop_defs.txt", col.names = FALSE,
            quote = FALSE, row.names = FALSE)

meta_df$id <- meta_df$id %>% 
  gsub("whtstbk_gbs_2012_brds_", "2012_", .) %>%
  gsub("\\.", "_", .) %>%
  str_trim

# genlight to genind (hack, genind interprets 0s as NAs instead of REF)
wht_gen_df <- as.data.frame(wht_gen) %>% data.table
rm(wht_gen)
wht_gen_df <- wht_gen_df + 1
wht_gen_df[is.na(wht_gen_df)] <- 0

# convert to genind
wht_genind <- df2genind(wht_gen_df, ncode = 1)
indNames(wht_genind) <- meta_df$id
save(wht_genind, file = "data/adegenet/wht_bi_snps_genind.R")
rm(wht_genind)

# convert to hierfstat
wht_fstat <- genind2hierfstat(wht_genind, pop = meta_df$pop)
save(wht_fstat, file = "data/adegenet/wht_bi_snps_hierfstat.R")

wht_basic <- basic.stats(wht_fstat)
save(wht_basic, file = "data/adegenet/hierfstat_basic.R")

# start here

load("data/adegenet/wht_bi_snps_hierfstat.R")
load("data/adegenet/hierfstat_basic.R")
load("data/adegenet/wht_bi_snps_genind.R")

wht_genpop <- genind2genpop(wht_genind, pop = meta_df$pop)
save(wht_genind, file = "data/adegenet/wht_bi_snps_genpop.R")
write.table(wht_genpop , file = "data/hierfstat/wht_genpop.gen")


