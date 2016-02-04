
# pacakge installation
# library(devtools)
# source("http://bioconductor.org/biocLite.R")
# biocLite("gdsfmt")
# biocLite("SNPRelate")
# install_github("slowkow/ggrepel")

# libraries
library(ggplot2)
library(ggrepel)
library(gdsfmt)
library(SNPRelate)
library(data.table)
library(dplyr)

list.files("functions", full.names = TRUE) %>% sapply(.,source, verbose = FALSE, echo = FALSE) %>% invisible


meta_df <- read.csv("metadata/mega_meta.csv")

# raw snps (raw ped file)
raw_snps <- data.frame(fread("data/other_formats/whtstbk_master_raw_no_sex.raw", stringsAsFactors = FALSE, header = TRUE))
#raw_snps <- read.table("data/other_formats/whtstbk_master_raw_no_sex.raw", stringsAsFactors = FALSE, header = TRUE)


geno_matrix <- data.matrix(raw_snps[,-c(1:6)])
sample_id <- raw_snps$FID
snp_id <- names(raw_snps[,-c(1:6)])

pos_df <- snp_id %>% fstat_label_to_columns()

snpgdsCreateGeno("whtstbk_raw_no_sex.gds", genmat = geno_matrix, 
                 sample.id = sample_id, snpfirstdim = FALSE, 
                 snp.id = snp_id, snp.chromosome = pos_df$chr, snp.position = pos_df$pos)

genofile <- snpgdsOpen("whtstbk_raw_no_sex.gds")

#snpgdsClose(genofile)

pca_samples <- sample_id[!grepl("2012", sample_id)]

pca <- snpgdsPCA(genofile, num.thread = 3, eigen.cnt = 6, missing.rate = 0.9, maf = 0.05, 
                 sample.id = pca_samples)

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
  geom_text(size = 4)



# extract the loadings for each snps
pos_df <- snp_id %>% fstat_label_to_columns()
pca_load <- snpgdsPCASNPLoading(pca, genofile)
pos_df$load <- pca_load$snploading[1,] 

pos_df <- pos_df %>% 
  mutate(pos = as.numeric(as.character(pos))) %>%
  arrange(chr, pos)

pos_df %>%
  ggplot(aes(x = pos, y = load))+
  geom_point()+
  geom_smooth()+
  facet_wrap(~chr)
