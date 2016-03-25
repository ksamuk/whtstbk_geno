
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

# Define palette

################################################################################
# plots
################################################################################

# IBS strip plot: morphotype comparison vs. year comparison

year_species <- diss_df %>%
  mutate(cluster_type = ifelse(cluster_1 == cluster_2, "Within Morphotype", "Between Morphotype")) %>%
  group_by(year_type, cluster_type, id) %>%
  summarise(mean_diss = mean(diss)) %>% ungroup %>%
  mutate(year_type = factor(year_type)) %>%
  mutate(cluster_type = factor(cluster_type)) %>%
  filter((year_type == "same" & cluster_type =="Between Morphotype") |(year_type == "different" & cluster_type =="Within Morphotype")) %>%
  ggplot(aes(x = cluster_type, y = mean_diss))+
  geom_jitter(color = "grey")+
  stat_summary(fun.y=median, fun.ymin=median, fun.ymax=median, 
               geom="crossbar", width=0.9, color = "black", fatten = 3)+
  theme_minimal() +
  theme(panel.grid.major.x = element_blank(),
        legend.position = "none",
        panel.margin = grid::unit(0, "cm"),
        axis.title.y=element_text(margin=margin(0,20,0,0), face = "bold"),
        axis.title.x=element_text(margin=margin(20,0,0,0), face = "bold"))+
  ylab("Genetic Similarity\n(% Shared SNPs)")+
  xlab("Comparison Type")+
  scale_x_discrete(breaks=c("Between Morphotype", "Within Morphotype"),
                   labels=c("Within Year\nWhite vs. Common", "Between Year\n White vs. White"))+
  ylim(0.74, 0.775)

sex_species <- diss_df %>%
  filter(year_1 == 2014 & year_2 ==2014) %>%
  #filter(cluster_1 == "wht" | cluster_2 == "wht") %>%
  #mutate(cluster_type = ifelse(cluster_1 == "wht" & cluster_2 == "wht", "White vs.\n White", ifelse(cluster_1 == "cmn" & cluster_2 == "cmn", "Common vs.\n Common", "White vs.\n Common"))) %>%
  mutate(cluster_type = ifelse(cluster_1 == cluster_2, "Within\n Morphotype", "Between\n Morphotype")) %>%
  group_by(cluster_type,sex_type, id) %>%
  summarise(mean_diss = mean(diss)) %>% ungroup %>%
  mutate(sex_type = factor(sex_type)) %>%
  mutate(cluster_type = factor(cluster_type)) %>%
  ggplot(aes(x = sex_type, y = mean_diss, color = sex_type))+
  #geom_dotplot(binaxis = "y", binwidth = 0.0001)+facet_grid(year_type~.)
  geom_jitter()+
  stat_summary(fun.y=median, fun.ymin=median, fun.ymax=median, 
               geom="crossbar", width=0.9, color = "black", fatten = 3)+
  #geom_point(position=position_jitterdodge(), alpha = 0.5) +
  theme_minimal() +
  facet_grid(.~cluster_type, switch = "x")+
  theme(panel.grid.major.x = element_blank(),
        #panel.grid.major.y = element_blank(),
        #legend.position = "bottom",
        #legend.direction = "vertical",
        legend.margin = grid::unit(0, "cm"),
        axis.text.x = element_blank(),
        panel.margin = grid::unit(0, "cm"),
        axis.title.y=element_text(margin=margin(0,20,0,0), face = "bold"),
        axis.title.x=element_text(margin=margin(20,0,0,0), face = "bold"))+
  scale_color_brewer("",labels = c("Between\nSexes\n", "Within\nSexes\n"), palette = "Set1") +
  ylab("\n")+
  xlab("Comparison Type")+
  ylim(0.74, 0.775)


fig_4 <- plot_grid(year_species, sex_species, align = "hv", labels = "AUTO", rel_widths = c(1, 1.2))
ggsave("figures/Figure4.pdf", plot = fig_4, width = 11, height = 8.5, scale = 0.8)

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
  #filter(cluster_1 == "wht" | cluster_2 == "wht") %>%
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

#result
# Coefficients:
# Estimate Std. Error  t value Pr(>|t|)    
#  (Intercept)                                 7.584e-01  1.680e-04 4513.698  < 2e-16 ***
#    cluster_typeWithin Morphotype               5.557e-03  2.376e-04   23.384  < 2e-16 ***
#    sex_typesame                                1.235e-03  2.376e-04    5.198 2.42e-07 ***
#    cluster_typeWithin Morphotype:sex_typesame -1.298e-05  3.361e-04   -0.039    0.969   
  
#                        Df    Sum Sq   Mean Sq   F value    Pr(>F)    
#  cluster_type             1 0.0081940 0.0081940 1091.0396 < 2.2e-16 ***
#  sex_type                 1 0.0004015 0.0004015   53.4640 5.182e-13 ***
#  cluster_type:sex_type    1 0.0000000 0.0000000    0.0015    0.9692 

################################################################################
# LN probe
################################################################################

diss_df %>%
  filter(cluster_1 != "cmn" & cluster_2 != "cmn") %>%
  filter(grepl("CL", id)&grepl("CL", id2)) %>% View
  mutate(is_ln_20 = (id == "CL50" | id2 == "CL50")) %>%
  ggplot(aes(x = cluster_type, y = diss, color = is_ln_20))+
  geom_boxplot()
