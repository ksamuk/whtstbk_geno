# per locus fst exploration

# libraries

#install.packages("qqman")
library("qqman")
library("dplyr")
library("ggplot2")
library("tidyr")
library("ggthemes")

list.files("functions", full.names = TRUE) %>% sapply(.,source, verbose = FALSE, echo = FALSE) %>% invisible

# read in raw data
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


#### READ IN XTX FILES FROM BAYPASS

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

#### Reading RSB / EHH files

rsb_df <- read.table("data/stats/whtstbk_rsb_df.txt", header = TRUE)
ies_df <- read.table("data/stats/whtstbk_ies_df.txt", header = TRUE)

####

#### Reading PCA loading file

pca_df <- read.table("data/stats/whtstbk_2014_pca_loadings.txt", header = TRUE)

####

# join all per locus files
fx_df <- full_join(fst_df, xtx_df) %>% full_join(rsb_df) %>% full_join(ies_df) %>% full_join(pca_df)
write.table(fx_df, "data/stats/snp_stats_master.txt", row.names = FALSE, quote = FALSE)

outlier_wht_cmn <- is.outlier(fx_df$fst_wht_cmn, cutoff = 0.99)
outlier_wht_cbr <- is.outlier(fx_df$fst_wht_cbr, cutoff = 0.99)
outlier_cbr_cmn <- is.outlier(fx_df$fst_cbr_cmn, cutoff = 0.99)

fx_df$outlier_fst_joint <- (outlier_wht_cmn | outlier_wht_cbr) & !(outlier_cbr_cmn)
fx_df$outlier_xtx_joint <- with(fx_df, (outlier_xtx_wht_cbr | outlier_xtx_wht_cmn) & !(outlier_xtx_cmn_cbr))

fst_long <- gather(fx_df, key = stat, value = value, -chr, -pos, -outlier_fst_joint, -outlier_xtx_joint)


fx_df %>%
  #sample_frac(0.01)%>%
  filter(chr == 18) %>%
  ggplot(aes(x = pos, y = ies_wht, color = "ies_wht")) +
  geom_point(size = 1)+
  geom_point(aes(y = -1*ies_cmn, color = "ies_cmn"), size = 1)+
  geom_smooth(aes(y = rsb_wht_cmn*100000, color = "rsb"))+
  #facet_grid(fst_type~chr, scales = "free", space = "free_x", switch = "x") +
  theme_classic()

fx_df %>%
  #sample_frac(0.01)%>%
  #filter(chr == 5) %>%
  ggplot(aes(x = rsb_wht_cmn, y = fst_wht_cmn)) +
  geom_point(size = 1)+
  #facet_grid(fst_type~chr, scales = "free", space = "free_x", switch = "x") +
  theme_classic()


fst_long %>%
  #filter(div_stat %in% c("wht_cmn", "wht_cbr", "cbr_cmn")) %>%
  #filter(chr %in% c(7)) %>%
  filter(!(grepl("outlier", stat))) %>% 
  filter(!(grepl("pval", stat))) %>% 
  filter((grepl("fst|rsb", stat))) %>%
  #sample_frac(0.01) %>%
  ggplot(aes(x = pos, y = value, color = outlier_xtx_joint)) +
  geom_point(size = 1)+
  facet_grid(stat~chr, scales = "free", space = "free_x", switch = "x") +
  theme_classic()+
  theme(axis.text.x = element_blank(), 
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        panel.margin = unit(0.1, "cm"),
        strip.background = element_blank())+
  scale_color_manual(values=c("black", "red"))

