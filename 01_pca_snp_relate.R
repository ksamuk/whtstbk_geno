
# pacakge installation
# library(devtools)
#source("http://bioconductor.org/biocLite.R")
#biocLite("gdsfmt")
#biocLite("SNPRelate")
# install_github("slowkow/ggrepel")

# libraries
library("ggplot2")
library("gdsfmt")
library("SNPRelate")
library("readr")
library("dplyr")
library("MASS")
library("fpc")

list.files("functions", full.names = TRUE) %>% sapply(.,source, verbose = FALSE, echo = FALSE) %>% invisible

select <- dplyr::select

meta_df <- read.csv("metadata/mega_meta.csv")

# raw snps (raw ped file)
#raw_snps <- data.frame(fread("data/snp_tables/whtstbk_all_pruned.gz", stringsAsFactors = FALSE, header = TRUE))
raw_snps <- read.table("data/snp_tables/whtstbk_2014_pruned.gz", header = TRUE, stringsAsFactors = FALSE)

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
      

snpgdsCreateGeno("data/snp_relate/whtstbk_2014_pruned.gds", genmat = geno_matrix, 
                 sample.id = sample_id, snpfirstdim = FALSE, 
                 snp.id = snp_id, snp.chromosome = pos_df$chr, snp.position = pos_df$pos)
# reclaim memory
rm(list=c("pos_df", "geno_matrix"))

# CAN START HERE
genofile_all <- snpgdsOpen("data/snp_relate/whtstbk_all_pruned.gds", allow.duplicate = TRUE)
genofile_2014 <- snpgdsOpen("data/snp_relate/whtstbk_2014_pruned.gds", allow.duplicate = TRUE)


#snpgdsClose(genofile)

pca_samples <- sample_id[!grepl("SK36|LN30|SR20|LN29|SR15", sample_id)]


diss <- snpgdsDiss(genofile_2014)
clust <- snpgdsHCluster(diss)
cut_tree <- snpgdsCutTree(clust)
tmp <- snpgdsDrawTree(clust, shadow.col = 10)

pca <- snpgdsPCA(genofile_2014, num.thread = 3, eigen.cnt = 16)


pca_load <- snpgdsPCASNPLoading(pca, genofile_2014, num.thread = 3)

load_df <- data.frame(snp = pca_load$snp.id, load_ev1 = t(pca_load$snploading)[,1], load_ev2 = t(pca_load$snploading)[,2])

load_df$chr <- load_df$snp %>% gsub("\\..*|chr", "",. ) %>% 
  gsub("Un", "XXII", .) %>% as.roman %>% as.numeric
  
load_df$pos <- load_df$snp %>% gsub("chr.*\\.|_[A-Z]*", "",. ) %>% as.numeric

load_df <- load_df %>%
  dplyr::select(chr, pos, load_ev1, load_ev2)

write.table(load_df, "data/stats/whtstbk_2014_pca_loadings.txt", row.names= FALSE, quote = FALSE)

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
#tab$id <- tab$id %>% gsub("whtstbk_gbs_|brds_", "", .)
pca_df <- left_join(tab, meta_df) %>%
  filter(!is.na(EV1))

pca_df$id <- pca_df$id %>% gsub("whtstbk_gbs_|brds_", "", .)
pca_df$id <- pca_df$id %>% gsub("NG-5241_[0-9]*_STD_", "", .)


pca_df %>%
  #filter(!(id %in% c("SR20", "LN29", "SR15", "2012_SF16"))) %>%
  filter(!is.na(cluster)) %>%
  #filter(region == "CB") %>%
  ggplot(aes(x = EV1, y = EV2, color = cluster, label = id))+
  geom_point(size = 2)
  #geom_text(size = 4)

pca_df %>%
  #filter(!(id %in% c("SR20", "LN29", "SR15", "2012_SF16"))) %>%
  filter(!is.na(cluster)) %>%
  #filter(region == "CB") %>%
  ggplot(aes(x = EV3, y = EV4, color = cluster, label = id))+
  geom_point(size = 2)
#geom_text(size = 4)

pca_df %>%
  filter(!(id %in% c("SR20", "LN29", "SR15"))) %>%
  #filter(region == "CB") %>%
  ggplot(aes(x = EV1, y = EV2, color = sex, label = id))+
  geom_point(size = 2)
#geom_text(size = 4)

parcoord(pca$eigenvect[,1:16], col = grepl("DK", pca_df$pop)+1, lty = 1)

parcoord(pca$eigenvect[,1:16], col = pca_df$sex, lty = 1)

kmeans_df <- pca_df[,2:10] %>%
  filter(!(is.na(EV1)))


# permutation of cluster separate

d <- dist(kmeans_df, method="euclidean") 

pfit <- hclust(d, method="ward.D")  

plot(pfit, labels = pca_df$id, cex = 0.2)   


# load the fpc package


# set the desired number of clusters                               
kbest.p <- 3      

#   Run clusterboot() with hclust 
#   ('clustermethod=hclustCBI') using Ward's method 
#   ('method="ward"') and kbest.p clusters 
#   ('k=kbest.p'). Return the results in an object 
#   called cboot.hclust.
cboot.hclust <- clusterboot(pmatrix, count = FALSE, B = 10000, 
                            clustermethod=kmeansCBI,
                            krange = 3, bootmethod=c("boot","noise","jitter"))

#Clusterwise Jaccard bootstrap (omitting multiple points) mean:
#  [1] 0.8923649 0.9042921 0.8957095

pca_df$sex<- as.factor(kmeans(kmeans_df , 2, iter.max = 1000, nstart = 100)$cluster)

# force MH5 into the GY cluster
cluster_labels <- c("M", "F")
pca_df$sex <- cluster_labels[match(pca_df$sex, c(1,2))]
#pca_df$cluster[pca_df$id =="MH5"] <- "cmn"

pca_df %>%
  #filter(region == "CB") %>%
  ggplot(aes(x = EV1, y = EV2, color = as.factor(sex), label = id))+
  geom_text(size = 5)

meta_df <- left_join(meta_df[,-5], pca_df[,c("id", "sex")])

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
