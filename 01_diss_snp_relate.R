
# pacakge installation
library(devtools)
#source("http://bioconductor.org/biocLite.R")
#biocLite("gdsfmt")
#biocLite("SNPRelate")
# install_github("slowkow/ggrepel")
install_github("dill/beyonce")

install_github("stanstrup/heatmap3")

################################################################################
# initials
################################################################################

# libraries
library("ggplot2")
library("gdsfmt")
library("SNPRelate")
library("readr")
library("dplyr")
library("tidyr")
library("MASS")
library("fpc")
library("ggthemes")
library("ggrepel")
library("RColorBrewer")
library("heatmap3")

list.files("functions", full.names = TRUE) %>% sapply(.,source, verbose = FALSE, echo = FALSE) %>% invisible

select <- dplyr::select

################################################################################
# raw data 
################################################################################

meta_df <- read.csv("metadata/mega_meta.csv")

# raw snps (raw ped file)
#raw_snps <- data.frame(fread("data/snp_tables/whtstbk_all_pruned.gz", stringsAsFactors = FALSE, header = TRUE))
raw_snps <- read.table("data/snp_tables/whtstbk_all_pruned.gz", header = TRUE, stringsAsFactors = FALSE)

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
pca_samples <- sample_id[!grepl("SK36|LN30|SR20|LN29|SR15|whtstbk_gbs_2012_brds_SR16", sample_id)]

snpgdsSampMissRate(genofile_all)

diss <- snpgdsIBS(genofile_all, sample.id = pca_samples, num.thread = 3, maf= 0.05, missing.rate = 0.05)
diss_mat <- diss$ibs
row.names(diss_mat) <- diss$sample.id
colnames(diss_mat) <- diss$sample.id

cols <- colorRampPalette(brewer.pal(10, "RdBu"))(256)

distCor <- function(x) as.dist(1-cor(t(x)))
hclustAvg <- function(x) hclust(x, method="average")

heatmap.3(diss_males, trace="none", zlim=c(-3,3), reorder=FALSE,
          distfun=distCor, hclustfun=hclustAvg, col=rev(cols), symbreak=FALSE) 


#clust <- snpgdsHCluster(diss)
#tmp <- snpgdsDrawTree(clust)

# males 2012 and 2014

male_ids <- meta_df %>%
  filter(sex == "M") %>%
  filter(pop %in% c("SR", "SF", "WR")) %>%
  select(id) %>% unlist %>% as.character

# remove low coverage individuals

diss_males <- diss_mat[colnames(diss_mat) %in% male_ids,row.names(diss_mat) %in% male_ids]

diss_males <- data.frame(diss_males)
diss_males$id <- row.names(diss_males)
diss_males <- gather(diss_males, key = id2, value = diss, -id)
diss_males <- left_join(diss_males, meta_df)


# create long version of matrix
diss_df_raw <- data.frame(diss_mat)
diss_df_raw$id <- row.names(diss_mat)

diss_df1 <- gather(diss_df_raw, key = id2, value = diss, -id)
names(diss_df1) <- c("id", "id2", "diss")
diss_df1 <- left_join(diss_df1, meta_df)

diss_df2 <- data.frame(diss_df1[,2])
names(diss_df2)[1] <- c("id")
diss_df2 <- left_join(diss_df2, meta_df)

diss_df <- data.frame(diss_df1[,c(1:3, 5, 8, 9)], diss_df2[, c(3,6:7)]) 
names(diss_df)[4:9] <- c("year_1", "cluster_1", "sex_1", "year_2", "cluster_2", "sex_2")

diss_df <- diss_df %>%
  mutate(year_type = ifelse(year_1 == year_2, "same", "different")) %>%
  mutate(cluster_type = ifelse(cluster_1 == cluster_2, "same", "different")) %>%
  mutate(sex_type = ifelse(sex_1 == sex_2, "same", "different"))

# Define palette

################################################################################
# plots
################################################################################

# heat map
myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")), space="Lab")

diss_df %>%
  ggplot(aes(x = id, y = id2, fill = diss))+
  geom_tile() +
  scale_fill_gradientn(colours = myPalette(100))+
  scale_y_discrete(expand = c(0, 0))+
  scale_x_discrete(expand = c(0, 0))+
  coord_equal()+
  theme_bw()+
  theme(axis.text=element_text(size=6),
        axis.title=element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1))

# strip plot

diss_df %>%
  group_by(year_type, cluster_type, id) %>%
  summarise(mean_diss = mean(diss)) %>% ungroup %>%
  mutate(year_type = factor(year_type)) %>%
  mutate(cluster_type = factor(cluster_type)) %>%
  ggplot(aes(x = year_type, fill = cluster_type, y = mean_diss, color = cluster_type))+
  geom_point(position=position_jitterdodge(), alpha = 0.5)+
  stat_summary(fun.y = 'mean', geom = 'errorbarh', 
               aes(xmin = (as.numeric(year_type)-as.numeric(cluster_type)), xmax = as.numeric(year_type)+as.numeric(cluster_type)), height = 0)

  ggplot(data = my.data, aes(x = tag2, y = values, color = as.factor(class))) +
  geom_jitter(position = position_jitter(width = .4)) +
   +
  scale_x_continuous(breaks = breaks, labels = levels(my.data$tag))
  
  diss_df %>%
    group_by(year_type, cluster_type, id) %>%
    summarise(mean_diss = mean(diss)) %>% ungroup %>%
    mutate(year_type = factor(year_type)) %>%
    mutate(cluster_type = factor(cluster_type)) %>%
  with(., as.numeric(year_type)-as.numeric(cluster_type)+1)





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
  geom_point(size = 3) +
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
