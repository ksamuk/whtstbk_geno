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
stats_folder <- "data/stats/2014"
stats_files <- list.files(stats_folder, pattern = ".txt", full.names = TRUE)

# names of the populations (cbr, wht, etc.)
pop_names <- c("cbr", "wht", "cmn")

################################################################################
# Summarise LD calculations
################################################################################

build_ld_decay_df <- function(pop_name, stats_files){
  
  # determine files
  stats_files_pop <- stats_files[grep(paste0(stats_folder,"/",pop_name), stats_files)]  
  stats_files_pop <- stats_files_pop[grep("r2", stats_files_pop)] 
  
  # read into a list df
  stats_dfs <- read.table(stats_files_pop[1], header = TRUE, stringsAsFactors = FALSE)
  
  # the magic: calculate site-wise ld (r2)
  # requires >3 pairwise r2's, each with >30 genotypes 
  
  
  sites_df <- stats_dfs %>%
    filter(N_INDV > 30) %>%
    filter(R.2 > 0.1) %>%
    group_by(CHR, POS1) %>%
    summarise(n_sites = n()) %>%
    filter(n_sites > 3) %>%
    select(CHR, POS1) %>%
    ungroup
  
  dist_df <- left_join(sites_df, stats_dfs, by = c("CHR", "POS1")) %>%
    mutate(dist = abs(POS1 - POS2)) %>%
    select(CHR, dist, R.2)
  
  names(dist_df) <- c("chr", "dist", paste0(pop_name, "_r2"))
  
  dist_df <- data.frame(pop = pop_name, dist_df[,c(1:3)])
  
  dist_df
    
}

wht_ld <- build_ld_decay_df("wht", stats_files)
cmn_ld <- build_ld_decay_df("cmn", stats_files)
cbr_ld <- build_ld_decay_df("cbr", stats_files)

ld_df <- full_join(wht_ld, cmn_ld)
ld_df <- full_join(ld_df, cbr_ld)

rm(wht_ld)
rm(cmn_ld)
rm(cbr_ld)

ld_df <- gather(ld_df, key = r2_type, value = r2, -pop, -chr, -dist)

ld_df <- ld_df %>%
  filter(!is.na(r2))%>%
  filter(!is.nan(r2))

ld_df$chr <- ld_df$chr %>% gsub("chr", "", .) %>% gsub("Un", "XXII", .) %>% as.roman %>% as.numeric

binomial_smooth <- function(...) {
  geom_smooth(method = "glm", method.args = list(family = "binomial"), ...)
}

ld_df %>%
  sample_frac(0.2)%>%
  #filter(r2_type == "wht_r2") %>%
  ggplot(aes(x = dist, y = r2, color = r2_type)) +
  stat_smooth(span = 0.01, method = "loess")+
  #binomial_smooth(se = FALSE)+
  facet_wrap(~chr)

