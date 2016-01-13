#library("adegenet")
library("hierfstat")
library("data.table")
library("dplyr")
library("stringr")

list.files("functions", full.names = TRUE) %>% sapply(.,source, verbose = FALSE, echo = FALSE) %>% invisible


# load hierfstat data
# 'wht_fstat'
load("data/adegenet/wht_bi_snps_hierfstat.R")

# read in metadata
meta_df <- read.csv("metadata/mega_meta.csv")
wht_fstat$pop

# calculate fst per locus
tmp <- fst_per_locus(wht_fstat[,1:1000], "AL", "BR")


