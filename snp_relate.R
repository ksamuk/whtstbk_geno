
# pacakge installation
# library(devtools)
#source("http://bioconductor.org/biocLite.R")
#biocLite("gdsfmt")
#biocLite("SNPRelate")
# install_github("slowkow/ggrepel")

# libraries
library(ggplot2)
library(gdsfmt)
library(SNPRelate)
library(data.table)
library(dplyr)
library(MASS)

list.files("functions", full.names = TRUE) %>% sapply(.,source, verbose = FALSE, echo = FALSE) %>% invisible

select <- dplyr::select

meta_df <- read.csv("metadata/mega_meta.csv")

# raw snps (raw ped file)
raw_snps <- data.frame(fread("data/other_formats/whtstbk_bial_no_sex.raw", stringsAsFactors = FALSE, header = TRUE))
#raw_snps <- read.table("data/other_formats/whtstbk_master_raw_no_sex.raw", stringsAsFactors = FALSE, header = TRUE)

# extract genotypes into matrix
geno_matrix <- data.matrix(raw_snps[,-c(1:6)])
sample_id <- raw_snps$FID
snp_id <- names(raw_snps[,-c(1:6)])

# reclaim memory
rm(raw_snps)

# fix formatting irregularities
pos_df <- snp_id %>% fstat_label_to_columns()
pos_df <- pos_df %>%
  mutate(chr = gsub("Un", "XXII", chr)) %>%
  mutate(chr = gsub("chr", "", chr) %>% as.character %>% as.roman %>% as.integer) %>%
  mutate(pos = pos %>% as.character %>% as.integer)
      

snpgdsCreateGeno("data/snp_relate/whtstbk_raw_no_sex_outgroup.gds", genmat = geno_matrix, 
                 sample.id = sample_id, snpfirstdim = FALSE, 
                 snp.id = snp_id, snp.chromosome = pos_df$chr, snp.position = pos_df$pos)
# reclaim memory
rm(list=c("pos_df", "geno_matrix"))

# CAN START HERE
genofile <- snpgdsOpen("data/snp_relate/whtstbk_raw_no_sex_outgroup.gds")

#snpgdsClose(genofile)

pca_samples <- sample_id[!grepl("SK36|2012|LN30|NG", sample_id)]

snpset <- snpgdsLDpruning(genofile, ld.threshold = 0.2, sample.id = pca_samples, 
                          missing.rate = 0.99, slide.max.n = 10, slide.max.bp = 10000)

snp_id <- snpset %>% unlist

snpgdsCreateGenoSet("data/snp_relate/whtstbk_raw_no_sex_outgroup.gds", "data/snp_relate/whtstbk_pruned.gds", 
                    snp.id = snp_id, sample.id = pca_samples)

pruned_set <- snpgdsOpen("data/snp_relate/whtstbk_pruned.gds", "data/other")
snpgdsGDS2PED(pruned_set, format = "1/2")

pca <- snpgdsPCA(genofile, num.thread = 3, eigen.cnt = 16, 
                 sample.id = pca_samples, snp.id = snp_id )

tab <- data.frame(id = pca$sample.id,
                  EV1 = pca$eigenvect[,1],    # the first eigenvector
                  EV2 = pca$eigenvect[,2],    # the second eigenvector
                  EV3 = pca$eigenvect[,3],    # the second eigenvector
                  EV4 = pca$eigenvect[,4],    # the second eigenvector
                  EV5 = pca$eigenvect[,5],    # the second eigenvector
                  EV6 = pca$eigenvect[,6],    # the second eigenvector
                  EV7 = pca$eigenvect[,7],    # the second eigenvector
                  stringsAsFactors = FALSE)
head(tab)

#convert ids to match metadata
tab$id <- tab$id %>% gsub("whtstbk_gbs_|brds_", "", .)
pca_df <- left_join(tab, meta_df) %>%
  filter(!is.na(EV1))

pca_df$id <- pca_df$id %>% gsub("NG-5241_[0-9]*_STD_", "", .)


pca_df %>%
  #filter(region == "CB") %>%
  ggplot(aes(x = EV1, y = EV2, color = region, label = id))+
  geom_text(size = 3)

parcoord(pca$eigenvect[,1:16], col = grepl("DK", pca_df$pop)+1, lty = 1)

parcoord(pca$eigenvect[,1:16], col = pca_df$year, lty = 1)

kmeans_df <- pca_df[,2:3] %>%
  filter(!(is.na(EV1)))

pca_df$cluster <- as.factor(kmeans(kmeans_df , 3, iter.max = 1000, nstart = 100)$cluster)

# force MH5 into the GY cluster
cluster_labels <- c("cmn", "cbr", "wht")
pca_df$cluster <- cluster_labels[match(pca_df$cluster, c(1,2,3))]
#pca_df$cluster[pca_df$id =="MH5"] <- "cmn"

pca_df %>%
  #filter(region == "CB") %>%
  ggplot(aes(x = EV1, y = EV2, color = as.factor(cluster), label = id))+
  geom_text(size = 5)

meta_df <- left_join(meta_df[,-8], pca_df[,c("id", "cluster")])

write.csv(meta_df, file = "metadata/mega_meta.csv", row.names = FALSE)

# extract the loadings for each snps

load_df <- snp_id %>% as.character %>% fstat_label_to_columns
pca_load <- snpgdsPCASNPLoading(pca, genofile)
load_df <- load_df %>%
  mutate(chr = gsub("Un", "XXII", chr)) %>%
  mutate(chr = gsub("chr", "", chr) %>% as.character %>% as.roman %>% as.integer)
  
load_df$load <- pca_load$snploading[1,] 

load_df <- load_df%>%
  mutate(pos = as.character(pos) %>% as.numeric) %>%
  arrange(chr, pos)

load_df %>%
  ggplot(aes(x = pos, y = abs(load)))+
  geom_point(size = 0.5)+
  #geom_smooth()+
  facet_wrap(~chr)

pca_df$region[83] <- "HA"

snpgdsFst(genofile, sample.id = , population = pca_df$region, method="W&H02")
