
# pacakge installation
#library(devtools)
#source("http://bioconductor.org/biocLite.R")
#biocLite("gdsfmt")
#biocLite("SNPRelate")
# install_github("slowkow/ggrepel")
#install_github("dill/beyonce")

#install_github("stanstrup/heatmap3")

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
# perform IBS analysis
################################################################################

# odd samples to remove (e.g. low data)
#LN30|LN29|
bad_samples <- c("SK36|SR20|SR15|whtstbk_gbs_2012_brds_SR16|whtstbk_gbs_2012_brds_SR12|whtstbk_gbs_2012_brds_PP8|whtstbk_gbs_2012_brds_CP21|CP23")
pca_samples <- sample_id[!grepl(bad_samples, sample_id)]

diss <- snpgdsIBS(genofile_all, sample.id = pca_samples, num.thread = 3, maf= 0.05, missing.rate = 0)
diss_mat <- diss$ibs
row.names(diss_mat) <- diss$sample.id
colnames(diss_mat) <- diss$sample.id

male_df <- meta_df %>%
  filter(sex == "M") %>%
  #filter(region !="CB") %>%
  filter(!(grepl(bad_samples, id))) 

male_ids <- male_df %>%  
  select(id) %>% 
  unlist %>% 
  as.character

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

# remove self comparisons
diss_df <- diss_df %>%
  filter(id != id2)

write.table(diss_df, "metadata/diss_df.txt", quote = FALSE, row.names = FALSE)

################################################################################
# stats
################################################################################

year_clust <- diss_df %>%
  #filter(cluster_1 == "wht" | cluster_2 == "wht") %>%
  #mutate(cluster_type = ifelse(cluster_1 == "wht" & cluster_2 == "wht", "White vs.\n White", ifelse(cluster_1 == "cmn" & cluster_2 == "cmn", "Common vs.\n Common", "White vs.\n Common"))) %>%
  mutate(cluster_type = ifelse(cluster_1 == cluster_2, "within_morpho", "between_morpho")) %>%
  group_by(year_type, cluster_type, id) %>%
  summarise(mean_diss = mean(diss)) %>% ungroup %>%
  filter((year_type == "same" & cluster_type =="between_morpho") |(year_type == "different" & cluster_type =="within_morpho")) %>%
  mutate(comparison = paste0(year_type,"_",cluster_type))


sex_clust <- diss_df %>%
  filter(year_1 == 2014 & year_2 ==2014) %>%
  filter(sex_1 == sex_2) %>%
  mutate(sex_type = ifelse(sex_1 == "F", "Female", "Male")) %>%
  #mutate(cluster_type = ifelse(cluster_1 == "wht" & cluster_2 == "wht", "White vs.\n White", ifelse(cluster_1 == "cmn" & cluster_2 == "cmn", "Common vs.\n Common", "White vs.\n Common"))) %>%
  mutate(cluster_type = ifelse(cluster_1 == cluster_2, "Within Morphotype", "Between Morphotype")) %>%
  group_by(cluster_type,sex_type, id) %>%
  summarise(mean_diss = mean(diss)) %>% ungroup

# linear model 1
year_clust %>%
  lm(data = ., mean_diss ~ comparison) %>%
  anova

# RESULT:
#Coefficients:
#  Estimate Std. Error t value Pr(>|t|)    
#(Intercept)                    0.7611513  0.0002008  3789.8   <2e-16 ***
#  comparisonsame_between_morpho -0.0028325  0.0002646   -10.7   <2e-16 ***
#comparison   1 0.0011971 0.00119710  114.58 < 2.2e-16 ***

# linear model 2
sex_clust %>%
  lm(data = ., mean_diss ~ cluster_type*sex_type) %>%
  anova()

# Coefficients:
#   Estimate Std. Error  t value Pr(>|t|)    
# (Intercept)                                7.597e-01  2.337e-04 3249.848   <2e-16 ***
#   cluster_typeWithin Morphotype              5.253e-03  3.306e-04   15.890   <2e-16 ***
#   sex_typeMale                               2.637e-05  3.318e-04    0.079    0.937    
# cluster_typeWithin Morphotype:sex_typeMale 5.863e-04  4.693e-04    1.249    0.212 

# Analysis of Variance Table
# 
# Response: mean_diss
# Df    Sum Sq   Mean Sq  F value Pr(>F)    
# cluster_type            1 0.0040874 0.0040874 558.2561 <2e-16 ***
#   sex_type                1 0.0000136 0.0000136   1.8544 0.1739    
# cluster_type:sex_type   1 0.0000114 0.0000114   1.5610 0.2121    
# Residuals             528 0.0038659 0.0000073                    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
