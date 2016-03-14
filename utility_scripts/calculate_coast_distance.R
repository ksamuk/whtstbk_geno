 ############################################################
# calculate euclidian & 'least cost' distance between populations
# KS Aug 2015
############################################################

rm(list =ls())

################################################################################
# Libraries and functions
################################################################################

library(parallel)
library(marmap)
library(dplyr)
library(fossil)

list.files("functions", full.names = TRUE) %>% sapply(source) %>% invisible

select <- dplyr::select

################################################################################
# Input files
################################################################################

pop.dat <- read.csv("metadata/whtstbk_site_coordinates.csv", header = TRUE, stringsAsFactors = FALSE)
fst.dat <- read.table("data/stats/fst_pairwise.txt", header = TRUE)

# shorten names for matching to geo distances file
fst.dat <- fst.dat %>%
  mutate(pop1_long = pop1) %>%
  mutate(pop2_long = pop2) %>%
  mutate(pop1 = gsub("_[a-z]*", "", pop1)) %>%
  mutate(pop2 = gsub("_[a-z]*", "", pop2)) %>%
  select(pop1, pop2, pop1_long, pop2_long, everything())

################################################################################
# Initialize marmap bathy data 
################################################################################

# Create nice looking color palettes
blues <- c("lightsteelblue4", "lightsteelblue3", "lightsteelblue2", "lightsteelblue1")
greys <- c(grey(0.4), grey(0.6), grey(0.8))

# preload and transform bathy data
bathy1<- getNOAA_bathy_prefix(-59.5, -65, 44, 47.5, res = 1, keep = TRUE, prefix = "metadata/")
trans.1.file <- "metadata/bathy.1.trans"

# if the transformed files don't exist, create them and save them to disk
if(!file.exists(trans.1.file)){
	trans1 <- trans.mat(bathy1, min.depth = 20, max.depth = -800)
	save(trans1, file = trans.1.file)
} else{
	load(trans.1.file )
}


################################################################################
# Run core function (compute_lc_distance_stats_file)
################################################################################


compute_lc_distance_stats_file(row_num, fst.dat, bathy1 = bathy1, trans1 = trans1, pop.dat = pop.dat)

row_nums <- 1:nrow(fst.dat)
distances.df <- lapply(row_nums, compute_lc_distance_stats_file, 
											 fst.dat = fst.dat, pop.dat = pop.dat, bathy1 = bathy1, trans1 = trans1, plot.map = FALSE)

# clean up results and write to file
distances.df <- bind_rows(distances.df)

distances.df <- distances.df %>%
	distinct

fst_dist <- cbind(fst.dat, distances.df[,-c(1:2)])

write.table(fst_dist, file = "metadata/fst_distances.txt", 
            quote = FALSE, row.names = FALSE)

fst_dist %>% 
  filter(type %in% c("wht_cmn", "wht_wht", "wht_cbr")) %>%
  mutate(tot_dist = least.cost.distance + dist.to.coast1 + dist.to.coast2) %>%
  #mutate(type = ifelse(type == "wht_wht", "wht_wht", "wht_cmn"))
  ggplot(aes(x = euc.distance, y = fst, color = type)) +
  geom_point()+
  geom_smooth(method = "lm", se = FALSE)

write.table(distances.df, file = "meta_data/pop_geo_distances.txt", 
						quote = FALSE, row.names = FALSE)
