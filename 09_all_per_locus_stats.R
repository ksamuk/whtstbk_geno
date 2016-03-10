# per locus fst exploration

################################################################################
# initials
################################################################################

library("qqman")
library("dplyr")
library("ggplot2")
library("tidyr")
library("ggthemes")

list.files("functions", full.names = TRUE) %>% sapply(.,source, verbose = FALSE, echo = FALSE) %>% invisible

# the raw "stats" files
stats_folder <- "data/stats/2014"
stats_files <- list.files(stats_folder, pattern = ".txt", full.names = TRUE)

################################################################################
# read in Fst data
################################################################################

perloc_files <- list.files("data/stats", pattern = "perloc", full.names = TRUE)

cbr_cmn_fst <- read.table(perloc_files[1], header = TRUE, stringsAsFactors = FALSE)
names(cbr_cmn_fst) <- c("chr", "pos", "fst")
cbr_cmn_fst$fst[cbr_cmn_fst$fst < 0] <- 0

wht_cmn_fst <- read.table(perloc_files[3], header = TRUE, stringsAsFactors = FALSE)
names(wht_cmn_fst) <- c("chr", "pos", "fst")
wht_cmn_fst$fst[wht_cmn_fst$fst < 0] <- 0


wht_cbr_fst <- read.table(perloc_files[2], header = TRUE, stringsAsFactors = FALSE)
names(wht_cbr_fst) <- c("chr", "pos", "fst")
wht_cbr_fst$fst[wht_cbr_fst$fst < 0] <- 0

fst_df <- full_join(wht_cmn_fst, wht_cbr_fst, by = c("chr", "pos"))
fst_df <- full_join(fst_df, cbr_cmn_fst, by = c("chr", "pos"))
names(fst_df)[3:5] <- c("fst_wht_cmn", "fst_wht_cbr", "fst_cbr_cmn")
fst_df$chr <- fst_df$chr %>% gsub("chr", "", .) %>% gsub("Un", "XXII", .) %>% as.roman %>% as.numeric

################################################################################
# READ IN XTX FILES FROM BAYPASS
################################################################################

slugs <- c("wht_cmn", "wht_cbr", "cmn_cbr", "wht_cmn.2012")

read_xtx_file <- function(slug){
  
  baypass_sites <- read.table("metadata/baypass_sites.txt", header = TRUE, stringsAsFactors = FALSE)
  
  # the baypass output files
  xtx_file <- list.files(paste0("data/baypass/", slug), pattern = paste0(slug, "_summary_pi_xtx"),full.names = TRUE)
  xtx_df <- read.table(xtx_file, header = TRUE, stringsAsFactors = FALSE)
  names(xtx_df) <- tolower(names(xtx_df)) 
  xtx_df <- cbind(baypass_sites, xtx_df[,-1])
  
  pod_xtx_file <- list.files(paste0("data/baypass/",slug), pattern = paste0(slug, "_pods_summary_pi_xtx"),full.names = TRUE)
  pod_xtx_df <- read.table(xtx_file, header = TRUE, stringsAsFactors = FALSE)
  names(pod_xtx_df) <- tolower(names(pod_xtx_df)) 
  pod_xtx_df <- cbind(baypass_sites, pod_xtx_df[,-1])
  
  #compute the 1% threshold
  pod_thresh <- quantile(pod_xtx_df$m_xtx, probs=0.99)
  
  xtx_df$outlier_xtx <- xtx_df$m_xtx >= pod_thresh
  
  names(xtx_df)[7] <- paste0("xtx_", slug)
  names(xtx_df)[9] <- paste0("outlier_xtx_", slug)
  xtx_df <- xtx_df[,c(1:2, 7, 9)]
  xtx_df
}

xtx_list <- lapply(slugs, read_xtx_file)

xtx_df <- Reduce(function(x, y) merge(x, y, all=TRUE), xtx_list) %>%
  arrange(chr, pos)

################################################################################
# Reading RSB / EHH / iHS files
################################################################################

rsb_df <- read.table("data/stats/whtstbk_rsb_df.txt", header = TRUE)
ies_df <- read.table("data/stats/whtstbk_ies_df.txt", header = TRUE)
ihs_df <- read.table("data/stats/whtstbk_ihs_df.txt", header = TRUE)

################################################################################
# Reading PCA loading file
################################################################################

pca_df <- read.table("data/stats/whtstbk_2014_pca_loadings.txt", header = TRUE)

################################################################################
# Summarise LD calculations
################################################################################

build_ld_df <- function(pop_name, stats_files, window = FALSE){
  
  # determine files
  stats_files_pop <- stats_files[grep(paste0(stats_folder,"/",pop_name), stats_files)]  
  stats_files_pop <- stats_files_pop[grep("r2", stats_files_pop)] 
  
  # read into a list df
  stats_dfs <- read.table(stats_files_pop[1], header = TRUE, stringsAsFactors = FALSE)
  
  # the magic: calculate site-wise ld (r2)
  # requires >3 pairwise r2's, each with >30 genotypes 
  calc_sitewise_ld <- . %>%
    filter(N_INDV > 30) %>%
    mutate(dist = abs(POS1 - POS2)) %>%
    group_by(CHR, POS1) %>%
    summarise(r2 = mean(R.2), n_sites = n()) %>%
    filter(n_sites > 3) %>%
    ungroup
  
  ld_df <- calc_sitewise_ld(stats_dfs)
  
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

wht_ld <- build_ld_df("wht", stats_files, window = FALSE)
cmn_ld <- build_ld_df("cmn", stats_files, window = FALSE)
cbr_ld <- build_ld_df("cbr", stats_files, window = FALSE)

ld_df <- left_join(wht_ld, cmn_ld)
ld_df <- left_join(ld_df, cbr_ld)

names(ld_df)[2] <- "pos"

ld_df$chr <- ld_df$chr %>% gsub("chr", "", .) %>% gsub("Un", "XXII", .) %>% as.roman %>% as.numeric

################################################################################
# Join all per locus files
################################################################################

fx_df <- full_join(fst_df, xtx_df) %>% full_join(rsb_df) %>% 
  full_join(ies_df) %>% full_join(ihs_df) %>% full_join(pca_df) %>% full_join(ld_df)
write.table(fx_df, "data/stats/snp_stats_master.txt", row.names = FALSE, quote = FALSE)

fx_df <- read.table("data/stats/snp_stats_master.txt", h = T)

outlier_wht_cmn <- is.outlier(fx_df$fst_wht_cmn, cutoff = 0.99)
outlier_wht_cbr <- is.outlier(fx_df$fst_wht_cbr, cutoff = 0.99)
outlier_cbr_cmn <- is.outlier(fx_df$fst_cbr_cmn, cutoff = 0.99)

fx_df$outlier_fst_joint <- (outlier_wht_cmn | outlier_wht_cbr) & !(outlier_cbr_cmn)
fx_df$outlier_xtx_joint <- with(fx_df, (outlier_xtx_wht_cbr | outlier_xtx_wht_cmn) & !(outlier_xtx_cmn_cbr))
fx_df$delta_ihs_wht_cmn <-fx_df$ihs_wht - fx_df$ihs_cmn
fx_df$delta_ihs_wht_cbr <-fx_df$ihs_wht - fx_df$ihs_cbr
fx_df$delta_ihs_cbr_cmn <-fx_df$ihs_cbr - fx_df$ihs_cmn

fx_df$delta_r2_wht_cmn <-fx_df$r2_wht - fx_df$r2_cmn
fx_df$delta_r2_wht_cbr <-fx_df$r2_wht - fx_df$r2_cbr
fx_df$delta_r2_cbr_cmn <-fx_df$r2_cbr - fx_df$r2_cmn

fst_long <- gather(fx_df, key = stat, value = value, -chr, -pos, -outlier_fst_joint, -outlier_xtx_joint)


fx_df %>%
  #sample_frac(0.01)%>%
  #filter(chr == 6) %>%
  ggplot(aes(x = pos, y = fst_wht_cmn)) +
  geom_point()+
  #geom_ribbon()+
  #geom_point(aes(y = -1*ihs_wht, color = "ihs_cmn"), size = 1)+
  #geom_line(aes(y = rsb_wht_cmn, color = "rsb"))+
  
  
  facet_grid(.~chr, scales = "free", space = "free_x", switch = "x") +
  theme_classic()+
  theme(axis.text.x = element_blank(), 
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        panel.margin = unit(0.1, "cm"),
        strip.background = element_blank())+
  scale_color_manual(values=c("black", "red"))

fst_long %>%
   sample_frac(0.01)%>%
  #filter(chr == 5) %>%
  filter(grepl("delta", stat)) %>%
  ggplot(aes(x = outlier_fst_joint, y = value)) +
  geom_jitter()+
  facet_wrap(~stat, scales = "free") +
  theme_classic()+
  theme(axis.text.x = element_blank(), 
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        panel.margin = unit(0.1, "cm"),
        strip.background = element_blank())


fst_long %>%
  #filter(div_stat %in% c("wht_cmn", "wht_cbr", "cbr_cmn")) %>%
  filter(chr %in% c(5)) %>%
  filter(!(grepl("outlier", stat))) %>% 
  filter(!(grepl("pval", stat))) %>% 
  filter((grepl("fst|delta", stat))) %>%
  #sample_frac(0.01) %>%
  ggplot(aes(x = pos, ymax = value, ymin = 0)) +
  #geom_point(size = 1)+
  geom_ribbon()+
  facet_grid(stat~., scales = "free", space = "free_x", switch = "x") +
  theme_classic()+
  theme(axis.text.x = element_blank(), 
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        panel.margin = unit(0.1, "cm"),
        strip.background = element_blank())+
  scale_color_manual(values=c("black", "red"))


fst_long %>%
  #filter(div_stat %in% c("wht_cmn", "wht_cbr", "cbr_cmn")) %>%
  ggplot(aes(x = outlier_xtx_joint, y = value, color = outlier_xtx_joint)) +
  geom_jitter(size = 1)+
  facet_wrap(~stat, scales = "free") +
  theme_classic()+
  theme(axis.text.x = element_blank(), 
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        panel.margin = unit(0.1, "cm"),
        strip.background = element_blank())+
  scale_color_manual(values=c("black", "red"))

