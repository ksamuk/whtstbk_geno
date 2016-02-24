# process population level summary statistics

################################################################################
# Libraries
################################################################################

library("dplyr")
library("ggplot2")
library("tidyr")
library("ggthemes")
library("viridis")

list.files("functions", full.names = TRUE) %>% sapply(.,source, verbose = FALSE, echo = FALSE) %>% invisible

################################################################################
# Raw data
################################################################################

# the raw "stats" files
stats_files <- list.files("data/stats", pattern = ".txt", full.names = TRUE)

# names of the populations (cbr, wht, etc.)
pop_names <- c("cbr", "wht", "cmn")

build_stats_df <- function(pop_name, stats_files){
  
  stats_files_pop <- stats_files[grep(pop_name, stats_files)] 
  stats_files_pop <- stats_files_pop[grep("chr", stats_files_pop, invert = TRUE)] 
  
  stats_dfs <- lapply(stats_files_pop, function(x) read.table(x, header = TRUE, stringsAsFactors = FALSE))
  names(stats_dfs) <- strsplit(stats_files_pop, "/") %>% 
    lapply(function(x)x[[3]]) %>% 
    unlist %>% 
    gsub(".txt", "", .)
  
  # normalize tajd windows (yes its a hack)
  stats_dfs[[grep("taj", names(stats_dfs))]]$BIN_START <- stats_dfs[[grep("taj", names(stats_dfs))]]$BIN_START + 1
  
  # drop the N_SNPs column for one of the stats 
  # will need to add another if n > 2 stats
  stats_dfs[[2]] <- stats_dfs[[2]] %>%
    select(-N_SNPS)
  
  # fix names
  stats_df <- left_join(stats_dfs[[1]], stats_dfs[[2]])
  names(stats_df) <- c("chr", "pos1", "pos2", "n_snps", "pi", "tajd")
  names(stats_df)[4:6] <- paste0(names(stats_df)[4:6], paste0("_", pop_name))
  
  stats_df[,paste0("pi_", pop_name)][is.nan(stats_df[,paste0("pi_", pop_name)])] <- NA
  stats_df[,paste0("tajd_", pop_name)][is.nan(stats_df[,paste0("tajd_", pop_name)])] <- NA
  stats_df
}

# build individual stats dfs
wht_df <- build_stats_df("wht", stats_files)
cmn_df <- build_stats_df("cmn", stats_files)
cbr_df <- build_stats_df("cbr", stats_files)

# join into a master df
stats_df <- left_join(wht_df, cmn_df)
stats_df <- left_join(stats_df, cbr_df)

# make a long version (for ggplottin')
stats_long <- gather(stats_df, key = stat, value = value, -chr, -pos1, -pos2)
stats_long$species <- strsplit(stats_long$stat, "_") %>% lapply(., function(x) x[[length(x)]]) %>% unlist
stats_long$stat <- gsub("_[a-z]{3}$", "", stats_long$stat)


################################################################################
# Summarise LD calculations
################################################################################

build_ld_df <- function(pop_name, stats_files, window = FALSE){
  
  # determine files
  stats_files_pop <- stats_files[grep(pop_name, stats_files)] 
  stats_files_pop <- stats_files_pop[grep("chr", stats_files_pop)] 
  
  # read into a list df
  stats_dfs <- lapply(stats_files_pop, function(x) read.table(x, header = TRUE, stringsAsFactors = FALSE))
  
  # rename list items to chromosomes
  names(stats_dfs) <- strsplit(stats_files_pop, "_") %>% 
    lapply(function(x)x[[2]]) %>% 
    unlist
  
  # the magic: calculate site-wise ld (r2)
  # requires >3 pairwise r2's, each with >30 genotypes 
  calc_sitewise_ld <- . %>%
    filter(N_INDV > 30) %>%
    mutate(dist = abs(POS1 - POS2)) %>%
    group_by(CHR, POS1) %>%
    summarise(r2 = mean(R.2), n_sites = n()) %>%
    filter(n_sites > 3) %>%
    ungroup
  
  ld_df <- lapply(stats_dfs, calc_sitewise_ld) %>% 
    rbind_all() 
  
  if (window == TRUE){
    ld_df$w_pos2 <- (((ld_df$POS1 / 50000) %>% floor) + 1)*50000
    ld_df$w_pos1 <- ld_df$w_pos2 - 49999
    
    ld_df <- ld_df %>%
      group_by(CHR, w_pos1, w_pos2) %>%
      summarise(r2 = mean(r2)) %>%
      ungroup
    
    names(ld_df)[4] <- paste0("r2_", pop_name)
    names(ld_df)[1] <- "chr"
    names(ld_df)[2] <- "pos1"
    names(ld_df)[3] <- "pos2"
    ld_df 
    
  } else{
    names(ld_df)[3] <- paste0("r2_", pop_name)
    names(ld_df)[1] <- "chr"
    names(ld_df)[2] <- "pos1"
    
    ld_df <- ld_df %>%
      select(-n_sites) %>% 
      ungroup()
    ld_df
  }

}

wht_ld <- build_ld_df("wht", stats_files, window = TRUE)
cmn_ld <- build_ld_df("cmn", stats_files, window = TRUE)
cbr_ld <- build_ld_df("cbr", stats_files, window = TRUE)

ld_df <- left_join(wht_ld, cmn_ld)
ld_df <- left_join(ld_df, cbr_ld)

################################################################################
# FST calculations
################################################################################

read_fst_window_file <- function(file_name){
  fst_df <- read.table(file_name, header = TRUE, stringsAsFactors = FALSE)
  fst_df <- fst_df[,c(1:3, 5)]
  comparison_name <- file_name %>% 
    gsub(".*/", "", .) %>% 
    strsplit(split = "_") %>% 
    unlist %>% .[c(1,3)] %>%
    paste( collapse = "_") %>% paste0("fst_",.)
  names(fst_df) <- c("chr", "pos1", "pos2", comparison_name)
  fst_df
}

wht_cbr_fst <- read_fst_window_file("data/stats/wht_vs_cbr_50k.windowed.weir.fst")
wht_cmn_fst <- read_fst_window_file("data/stats/wht_vs_cmn_50k.windowed.weir.fst")
cbr_cmn_fst <- read_fst_window_file("data/stats/cbr_vs_cmn_50k.windowed.weir.fst")

fst_df <- left_join(wht_cbr_fst, wht_cmn_fst)
fst_df <- left_join(fst_df, cbr_cmn_fst)

################################################################################
# Join all data
################################################################################

stat_df <- left_join(ld_df, stats_df)
stat_df <- left_join(stat_df, fst_df)

# joint outlier score

wht_cbr_outlier <- stat_df$fst_wht_cbr %>% is.outlier
wht_cmn_outlier <- stat_df$fst_wht_cmn %>% is.outlier
cbr_cmn_outlier <- stat_df$fst_cbr_cmn %>% is.outlier

stat_df$fst_joint_outlier <- (wht_cbr_outlier & wht_cmn_outlier) & !(cbr_cmn_outlier)

write.table(stat_df, "data/pop_gen_stats_50k.txt", quote = FALSE, row.names = FALSE)
stat_df <- read.table("data/pop_gen_stats_50k.txt", header = TRUE, stringsAsFactors = FALSE)

stat_long <- gather(stat_df, key = stat, value = value, -chr, -pos1, -pos2)

################################################################################
# Plots
################################################################################

# > names(stat_df)
# [1] "chr"         "pos1"        "pos2"        "r2_wht"      "r2_cmn"      "r2_cbr"      "n_snps_wht" 
# [8] "pi_wht"      "tajd_wht"    "n_snps_cmn"  "pi_cmn"      "tajd_cmn"    "n_snps_cbr"  "pi_cbr"     
# [15] "tajd_cbr"    "wht_cbr_fst" "wht_cmn_fst" "cbr_cmn_fst"

stat_df %>%
  filter(chr == "chrI") %>%
  ggplot(aes(x = pos1, y = wht_cmn_fst %>% scale, color = "fst"))+
  geom_smooth(se = FALSE)+
  geom_smooth(aes(y = wht_cbr_fst %>% scale, color = "fst"), se = FALSE)+
  geom_smooth(aes(y = r2_wht %>% scale, color = "wht_r2"), se = FALSE)+
  geom_smooth(aes(y = tajd_wht %>% scale, color = "wht_tajd"), se = FALSE)+
  geom_smooth(aes(y = pi_wht %>% scale, color = "wht_pi"), se = FALSE)+
  facet_grid(chr~.)


stat_df %>%
  ggplot(aes(x = pos1, y = tajd_wht, color = fst_joint_outlier))+
  geom_point()+
  facet_grid(chr~.)

  
stat_df %>%
  filter(!is.na(fst_joint_outlier)) %>%
  ggplot(aes(y = asin(sqrt(r2_wht)), x = chr, color = fst_joint_outlier))+
  geom_boxplot()+
  facet_grid(chr~.)
  scale_color_viridis()

stat_long %>%
  filter((grepl("tajd", stat))) %>%
  ggplot(aes(x = value))+
  geom_histogram()+
  #stat_smooth()+
  facet_grid(~stat, scales = "free_y")

stats_df %>%
  mutate(tajd_diff = abs(scale(tajd_wht, center = FALSE) - scale(tajd_cmn, center = FALSE))) %>%
  ggplot(aes(y = tajd_diff, x = tajd_wht, color = tajd_wht)) +
  geom_point()+
  geom_smooth(method = "lm")+
  tehe
  facet_wrap(~chr)

stats_df %>%
  ggplot(aes(y = tajd_wht, x = pos2)) +
  geom_smooth()+
  geom_smooth(aes(y = tajd_cmn))+
  geom_smooth(aes(y = tajd_cbr))+
  facet_wrap(~chr)

stats_df %>%
  ggplot(aes(y = pi_wht, x = pos2)) +
  geom_smooth()+
  geom_smooth(aes(y = pi_cmn))+
  geom_smooth(aes(y = pi_cbr))+
  facet_wrap(~chr)

stats_long %>%
  ggplot(aes(x = value)) +
  geom_histogram()+
  facet_wrap(species~stat, scales = "free")

stat_long %>%
  filter(chr == "chrVII") %>%
  #filter(stat == "fst_joint_outlier") %>%
  ggplot(aes(x = pos2, y = value)) +
  geom_line()+
  facet_grid(stat~., scales = "free_y")




