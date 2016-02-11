# prepare a LD-pruned snp file for fastStructure
# the key problem this solves is splitting a multipopulation snp file 
# into multiple files, applying the LD pruning, and reforming the multipop file

# Kieran Samuk Feb 10/2016

################################################################################
# Libraries
################################################################################

library(ggplot2)
library(gdsfmt)
library(SNPRelate)
library(data.table)
library(dplyr)
library(MASS)

list.files("functions", full.names = TRUE) %>% sapply(.,source, verbose = FALSE, echo = FALSE) %>% invisible
select <- dplyr::select

################################################################################
# Read in raw data
################################################################################

# the meta data file (pop codes)
meta_df <- read.csv("metadata/mega_meta.csv")

# raw snps (raw ped file)
raw_snps <- data.frame(fread("data/other_formats/whtstbk_bial_no_sex.raw", stringsAsFactors = FALSE, header = TRUE))

# remove 2012 data 
raw_snps <- raw_snps %>%
  filter(!grepl("2012", raw_snps$FID))

# extract genotypes into matrix
geno_matrix <- data.matrix(raw_snps[,-c(1:6)])

# grab names and convery to matching meta_df coding
sample_df <- data.frame(id = raw_snps$FID %>% gsub("whtstbk_gbs_|_brds", "", .))
sample_df <- left_join(sample_df, meta_df)
sample_df$row <- 1:length(sample_df$id) 

# grab snp names (for chr pos)
snp_id <- names(raw_snps[,-c(1:6)])

# reclaim memory
rm(raw_snps)

# formate chr pos data
pos_df <- snp_id %>% fstat_label_to_columns()
pos_df <- pos_df %>%
  mutate(chr = gsub("Un", "XXII", chr)) %>%
  mutate(chr = gsub("chr", "", chr) %>% as.character %>% as.roman %>% as.integer) %>%
  mutate(pos = pos %>% as.character %>% as.integer)


################################################################################
# Create a pruned .gds file for each region
################################################################################

# split by region
reg_df <- split(sample_df, sample_df$region)

lapply(reg_df, function(x) split(x, x$cluster))

# repair DK missing clusters
reg_df$DK$cluster <- reg_df$DK$species


snpgdsCreateGeno("data/snp_relate/pop/whtstbk_raw_no_sex_outgroup.gds", genmat = geno_matrix, 
                 sample.id = sample_id, snpfirstdim = FALSE, 
                 snp.id = snp_id, snp.chromosome = pos_df$chr, snp.position = pos_df$pos)

# CAN START HERE
genofile <- snpgdsOpen("data/snp_relate/tmp.gds")

# reclaim memory
rm(list=c("pos_df", "geno_matrix"))

#snpgdsClose(genofile)

pca_samples <- sample_id[!grepl("SK36|2012|LN30|NG", sample_id)]

snpset <- snpgdsLDpruning(genofile, ld.threshold = 0.2, sample.id = pca_samples, 
                          missing.rate = 0.99, slide.max.n = 10, slide.max.bp = 10000)