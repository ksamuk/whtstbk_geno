
# pacakge installation
# library(devtools)
#source("http://bioconductor.org/biocLite.R")
#biocLite("gdsfmt")
#biocLite("SNPRelate")
# install_github("slowkow/ggrepel")
# install_github("dill/beyonce")

################################################################################
# initials
################################################################################

# libraries
library("ggplot2")
library("gdsfmt")
library("SNPRelate")
library("readr")
library("dplyr")
library("MASS")
library("fpc")
library("ggthemes")
library("ggrepel")
library("RColorBrewer")
library("bigmemory")

list.files("functions", full.names = TRUE) %>% sapply(.,source, verbose = FALSE, echo = FALSE) %>% invisible

select <- dplyr::select

################################################################################
# raw data 
################################################################################

meta_df <- read.csv("metadata/mega_meta.csv")

# raw snps (raw ped file)
# read in geno matrix as big matrix
geno_matrix <- read.big.matrix("data/snp_tables/whtstbk_bial_nomaf_nosex.raw", sep=" ", header = TRUE)
geno_matrix  <- geno_matrix[,-c(1:6)]

# read in FIDs

fid_df <- read.table(pipe("cut -d \" \" -f1 data/snp_tables/whtstbk_bial_nomaf_nosex.raw"), header = TRUE)

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
      
snpgdsCreateGeno("data/snp_relate/whtstbk_all_nomaf_unpruned.gds", genmat = geno_matrix,
                sample.id = sample_id, snpfirstdim = FALSE,
                snp.id = snp_id, snp.chromosome = pos_df$chr, snp.position = pos_df$pos)

rm(list=c("pos_df", "geno_matrix"))

# CAN START HERE
genofile_all <- snpgdsOpen("data/snp_relate/whtstbk_all_nomaf_unpruned.gds", allow.duplicate = TRUE)
#genofile_2014 <- snpgdsOpen("data/snp_relate/whtstbk_2014_pruned.gds", allow.duplicate = TRUE)

#snpgdsClose(genofile)

################################################################################
# perform pca
################################################################################

# odd samples to remove (e.g. low data)
pca_samples <- sample_id[!grepl("SK36|LN30|SR20|LN29|SR15|SR131|whtstbk_gbs_2012_SR16", sample_id)]

#diss <- snpgdsDiss(genofile_all, sample.id = pca_samples)
#clust <- snpgdsHCluster(diss)
#tmp <- snpgdsDrawTree(clust)

# do the actual PCA
pca <- snpgdsPCA(genofile_all, num.thread = 3, eigen.cnt = 16, sample.id = pca_samples, maf=0.05, missing.rate=0.8,)

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

pca_df %>%
  mutate(Year = as.factor(year) %>% reorder (., . == "2012")) %>%
  mutate(Sex = sex) %>%
  mutate(Sex = ifelse(Sex == "M", "Male", "Female")) %>%
  mutate(Group = cluster_names[match(cluster, cluster_short)]) %>%
  mutate(Region = region_names[match(region, region_short)]) %>%
  #filter(!(id %in% c("SR20", "LN29", "SR15", "2012_SF16"))) %>%
  filter(!is.na(cluster)) %>%
  #filter(region == "CB") %>%
  ggplot(aes(x = EV1, y = EV2, color = Sex, label = id, shape = Year))+
  #geom_point(size = 3) +
  geom_text()+
  theme_bw() +
  theme(legend.key = element_blank())+
  xlab(paste0("PC1 ", "(",var_prop[1], "%)"))+
  ylab(paste0("PC2 ", "(",var_prop[2], "%)"))+
  #facet_wrap(~region, scales = "free")+
  #scale_color_manual(values = pal)
  scale_color_brewer(palette = "Set1")
 

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
