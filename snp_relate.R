
# pacakge installation
# library(devtools)
source("http://bioconductor.org/biocLite.R")
biocLite("gdsfmt")
biocLite("SNPRelate")
# install_github("slowkow/ggrepel")

# libraries
library(ggplot2)
library(gdsfmt)
library(SNPRelate)
library(data.table)
library(dplyr)
library(MASS)

list.files("functions", full.names = TRUE) %>% sapply(.,source, verbose = FALSE, echo = FALSE) %>% invisible


meta_df <- read.csv("metadata/mega_meta.csv")

# raw snps (raw ped file)
raw_snps <- data.frame(fread("data/other_formats/whtstbk_master_raw_no_sex.raw", stringsAsFactors = FALSE, header = TRUE))
#raw_snps <- read.table("data/other_formats/whtstbk_master_raw_no_sex.raw", stringsAsFactors = FALSE, header = TRUE)


geno_matrix <- data.matrix(raw_snps[,-c(1:6)])
sample_id <- raw_snps$FID
snp_id <- names(raw_snps[,-c(1:6)])

pos_df <- snp_id %>% fstat_label_to_columns()
pos_df <- pos_df %>%
  mutate(chr = gsub("Un", "XXII", chr)) %>%
  mutate(chr = gsub("chr", "", chr) %>% as.character %>% as.roman %>% as.integer) %>%
  mutate(pos = pos %>% as.character %>% as.integer)
      

snpgdsCreateGeno("whtstbk_raw_no_sex.gds", genmat = geno_matrix, 
                 sample.id = sample_id, snpfirstdim = FALSE, 
                 snp.id = snp_id, snp.chromosome = pos_df$chr, snp.position = pos_df$pos)

genofile <- snpgdsOpen("whtstbk_raw_no_sex.gds")

#snpgdsClose(genofile)

pca_samples <- sample_id[!grepl("2012|SK36", sample_id)]

snpset <- snpgdsLDpruning(genofile, ld.threshold = 0.2, sample.id = pca_samples, missing.rate = 0.05, maf = 0.05)

snp_id <- snpset %>% unlist

pca <- snpgdsPCA(genofile, num.thread = 3, eigen.cnt = 16, bayesian = TRUE, sample.id= pca_samples, snp.id = snp_id)

tab <- data.frame(id = pca$sample.id,
                  EV1 = pca$eigenvect[,1],    # the first eigenvector
                  EV2 = pca$eigenvect[,2],    # the second eigenvector
                  EV3 = pca$eigenvect[,3],    # the second eigenvector
                  EV4 = pca$eigenvect[,4],    # the second eigenvector
                  stringsAsFactors = FALSE)
head(tab)

#convert ids to match metadata
tab$id <- tab$id %>% gsub("whtstbk_gbs_|brds_", "", .)
pca_df <- left_join(meta_df, tab)


pca_df %>%
  ggplot(aes(x = EV1, y = EV2, color = region, label = species))+
  geom_text(size = 5)

parcoord(pca$eigenvect[,1:16], col = pca_df$region, lty = 1)


# extract the loadings for each snps

load_df <- snp_id %>% as.character %>% fstat_label_to_columns
pca_load <- snpgdsPCASNPLoading(pca, genofile)
load_df <- load_df %>%
  mutate(chr = gsub("Un", "XXII", chr)) %>%
  mutate(chr = gsub("chr", "", chr) %>% as.character %>% as.roman %>% as.integer) %>%
  
load_df$load <- pca_load$snploading[1,] 

load_df <- load_df%>%
  mutate(pos = as.character(pos) %>% as.numeric) %>%
  arrange(chr, pos)

load_df %>%
  ggplot(aes(x = pos, y = abs(load)))+
  geom_point(size = 0.5)+
  #geom_smooth()+
  facet_wrap(~chr)
