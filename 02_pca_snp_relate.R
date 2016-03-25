# pacakge installation
#library(devtools)
#source("http://bioconductor.org/biocLite.R")
#biocLite("gdsfmt")
#biocLite("SNPRelate")
# install_github("slowkow/ggrepel")
#install_github("dill/beyonce")

#install_github("leeper/colourlovers")

################################################################################
# initials
################################################################################

# libraries
library("ggplot2")
library("gdsfmt")
library("SNPRelate")
library("dplyr")
library("tidyr")
library("ggthemes")
library("RColorBrewer")
library("bigmemory")
library("cowplot")
library("viridis")
library("colourlovers")

list.files("functions", full.names = TRUE) %>% sapply(.,source, verbose = FALSE, echo = FALSE) %>% invisible

select <- dplyr::select

################################################################################
# raw data 
################################################################################

meta_df <- read.csv("metadata/mega_meta.csv")

# raw snps (raw ped file)
#raw_snps <- data.frame(fread("data/snp_tables/whtstbk_all_pruned.gz", stringsAsFactors = FALSE, header = TRUE))
geno_matrix <- read.big.matrix("data/snp_tables/whtstbk_all_pruned.raw", sep=" ", header = TRUE)
geno_matrix  <- geno_matrix[,-c(1:6)]

# read in FIDs

fid_df <- read.table(pipe("cut -d \" \" -f1 data/snp_tables/whtstbk_all_pruned.raw"), header = TRUE)

# extract genotypes into matrix
sample_id <- fid_df$FID
snp_id <- colnames(geno_matrix)

# fix formatting irregularities
fstat_label_to_columns <- function(x){
  if(!is.null(names(x))){
    x <- names(x)
  }
  chr <- gsub("\\:.*", "", x)
  pos <- gsub("[^0-9]", "", x)
  data.frame(chr, pos)
}
pos_df <- snp_id %>% fstat_label_to_columns()
pos_df <- pos_df %>%
  mutate(chr = gsub("Un", "XXII", chr)) %>%
  mutate(chr = gsub("chr", "", chr) %>% as.character %>% as.roman %>% as.integer) %>%
  mutate(pos = pos %>% as.character %>% as.integer)

################################################################################
# create/open gen object 
################################################################################

#snpgdsCreateGeno("data/snp_relate/whtstbk_all_pruned.gds", genmat = geno_matrix, 
#                 sample.id = sample_id, snpfirstdim = FALSE, 
#                 snp.id = snp_id, snp.chromosome = pos_df$chr, snp.position = pos_df$pos)
# reclaim memory
rm(list=c("pos_df", "geno_matrix"))

# CAN START HERE
genofile_all <- snpgdsOpen("data/snp_relate/whtstbk_all_pruned.gds", allow.duplicate = TRUE)
#genofile_2014 <- snpgdsOpen("data/snp_relate/whtstbk_2014_pruned.gds", allow.duplicate = TRUE)

#snpgdsClose(genofile)

################################################################################
# perform pca
################################################################################

# odd samples to remove (e.g. low data)
bad_samples <- c("SK36|LN30|SR20|LN29|SR15|whtstbk_gbs_2012_brds_SR16|whtstbk_gbs_2012_brds_SR12|whtstbk_gbs_2012_brds_PP8|whtstbk_gbs_2012_brds_CP21|CP23")
pca_samples <- sample_id[!grepl(bad_samples, sample_id)]


#diss <- snpgdsDiss(genofile_all, sample.id = pca_samples)
#clust <- snpgdsHCluster(diss)
#tmp <- snpgdsDrawTree(clust)

# do the actual PCA
pca <- snpgdsPCA(genofile_all, num.thread = 3, eigen.cnt = 16, sample.id = pca_samples, missing.rate = 0.01, maf = 0.05)
#pca <- snpgdsPCA(genofile_all, num.thread = 3, eigen.cnt = 16, missing.rate = 0.05, maf = 0.05)

# extract loadings
pca_load <- snpgdsPCASNPLoading(pca, genofile_all, num.thread = 3)
load_df <- data.frame(snp = pca_load$snp.id, load_ev1 = t(pca_load$snploading)[,1], load_ev2 = t(pca_load$snploading)[,2])
load_df$chr <- load_df$snp %>% gsub("\\..*|chr", "",. ) %>% 
  gsub("Un", "XXII", .) %>% as.roman %>% as.numeric
  
load_df$pos <- load_df$snp %>% gsub("chr.*\\.|_[A-Z]*", "",. ) %>% as.numeric
load_df <- load_df %>%
  dplyr::select(chr, pos, load_ev1, load_ev2)
write.table(load_df, "data/stats/whtstbk_2014_pca_loadings.txt", row.names= FALSE, quote = FALSE)

# create PCA data frame
tab <- data.frame(id = pca$sample.id,
                  EV1 = pca$eigenvect[,1],    # the first eigenvector
                  EV2 = pca$eigenvect[,2],    # the second eigenvector
                  EV3 = pca$eigenvect[,3],    # the second eigenvector
                  EV4 = pca$eigenvect[,4],    # the second eigenvector
                  EV5 = pca$eigenvect[,5],    # the second eigenvector
                  EV6 = pca$eigenvect[,6],    # the second eigenvector
                  EV7 = pca$eigenvect[,7],    # the second eigenvector
                  stringsAsFactors = FALSE)

#convert ids to match metadata
#tab$id <- tab$id %>% gsub("whtstbk_gbs_|brds_", "", .)
pca_df <- left_join(tab, meta_df) %>%
  filter(!is.na(EV1))

pca_df$id <- pca_df$id %>% gsub("whtstbk_gbs_|brds_", "", .)
pca_df$id <- pca_df$id %>% gsub("NG-5241_[0-9]*_STD_", "", .)

# percent variance explained by PCs
var_prop <- pca$varprop[1:2]
var_prop <- pca$varprop[1:2]*100
var_prop <- signif(var_prop, 2)

################################################################################
# plots
################################################################################

# color palatte

pal <- c(rev(brewer.pal(3, "Set1")), "#984EA3")
#pal <- rev(wes_palette("Royal2", 3, type = "continuous"))
#pal <- c("#4bbfc0","#e87348","#3f95ca","#c43896")

# clusters

# region names

cluster_short <- c("cmn", "cbr", "wht")
cluster_names <- c("Mainland Common", "Bras d'Or Common", "White")
region_short <- c("AN", "CB", "HA", "GY")
region_names <- c("Antigonish", "Bras d'Or", "Halifax", "Guysborough")
pop_names <- c("Antigonish Land", "Black River", "Canal Lake", "Captain's Pond", "Gillies Cove", "Little Narrows",
               "Milford Haven", "Middle River", "River Tillard", "St. Francais Harbour", "Sheet Harbour", "Skye River",
               "Salmon River", "Porper Pond", "Pomquet", "Right's River")

palette4 <- clpalette('694737')$colors %>% as.character %>% paste0("#", .)

figure_2 <- pca_df %>%
  mutate(Year = as.factor(year) %>% reorder (., . == "2012")) %>%
  mutate(year_pch = c(21,24)[as.numeric(Year)]) %>%
  mutate(Sex = sex) %>%
  mutate(Sex = ifelse(Sex == "M", "Male", "Female")) %>%
  mutate(Group = cluster_names[match(cluster, cluster_short)]) %>%
  mutate(Region = region_names[match(region, region_short)]) %>%
  #filter(!(id %in% c("SR20", "LN29", "SR15", "2012_SF16"))) %>%
  filter(!is.na(cluster)) %>%
  #filter(region == "CB") %>%
  ggplot(aes(x = EV1, y = EV2, fill = Sex, shape=Year))+
  geom_point(aes(fill=Sex, shape=Year), colour="black", size=5)+
  xlab(paste0("PC1 ", "(",var_prop[1], "%)"))+
  ylab(paste0("PC2 ", "(",var_prop[2], "%)"))+
  theme_bw() +
  theme(legend.key = element_blank(),
        text = element_text(size = 18),
        axis.title.y = element_text(margin=margin(0,20,0,0), face = "bold"),
        axis.title.x = element_text(margin=margin(20,0,0,0), face = "bold"))+
  scale_fill_fivethirtyeight()+
  scale_shape_manual(values = c(21,24))

ggsave(plot = figure_2, "figures/Figure2.pdf", width = 11, height = 8.5)

figure_S2 <- pca_df %>%
  mutate(Year = as.factor(year) %>% reorder (., . == "2012")) %>%
  mutate(year_pch = c(21,24)[as.numeric(Year)]) %>%
  mutate(Sex = sex) %>%
  mutate(Sex = ifelse(Sex == "M", "Male", "Female")) %>%
  mutate(Group = cluster_names[match(cluster, cluster_short)]) %>%
  mutate(Region = region_names[match(region, region_short)]) %>%
  #filter(!(id %in% c("SR20", "LN29", "SR15", "2012_SF16"))) %>%
  filter(!is.na(cluster)) %>%
  #filter(region == "CB") %>%
  ggplot(aes(x = EV1, y = EV2, fill = Group))+
  geom_point(aes(fill = Group), colour = "black", size=5, pch = 21)+
  xlab(paste0("PC1 ", "(",var_prop[1], "%)"))+
  ylab(paste0("PC2 ", "(",var_prop[2], "%)"))+
  theme_bw() +
  theme(legend.key = element_blank(),
        text = element_text(size = 18),
        axis.title.y = element_text(margin=margin(0,20,0,0), face = "bold"),
        axis.title.x = element_text(margin=margin(20,0,0,0), face = "bold"))+
  scale_fill_manual(values = c("#77AB43", "#008FD5", "#FFFFFF"))

ggsave(plot = figure_S2, "figures/FigureS2.pdf", width = 11, height = 8.5)
 

################################################################################
# test of cluster existance
################################################################################

# subset of the pcs for  cluster permutation 
kmeans_df <- pca_df[,2:3] %>%
  filter(!(is.na(EV1)))

# quick k-means

clusters <- kmeans(kmeans_df, 3, iter.max = 1000, nstart = 100)$cluster

cluster_df <- data.frame(pca_df, clusters)


# calculate disance matrix using pc data
d <- dist(kmeans_df, method="euclidean") 

# load the fpc package
library("fpc")

# run clusterboot on pca data
cboot.hclust <- clusterboot(d, count = TRUE, B = 1000, 
                            clustermethod=kmeansCBI,
                            krange = 3, bootmethod=c("boot","noise","jitter"))

# * Cluster stability assessment *
#   Cluster method:  kmeans 
# Full clustering results are given as parameter result
# of the clusterboot object, which also provides further statistics
# of the resampling results.
# Number of resampling runs:  1000 
# 
# Number of clusters found in data:  3 
# 
# Clusterwise Jaccard bootstrap (omitting multiple points) mean:
#   [1] 0.9203614 0.9066442 0.8857714
# dissolved:
#   [1]   0  48 165
# recovered:
#   [1] 807 790 794
# Clusterwise Jaccard replacement by noise mean:
#   [1] 0.9435032 0.9280240 0.9104470
# dissolved:
#   [1]   0  49 155
# recovered:
#   [1] 865 850 844
# Clusterwise Jaccard jittering mean:
#   [1] 0.9365534 0.9274288 0.9123906
# dissolved:
#   [1]   0  48 113
# recovered:
#   [1] 839 839 839
