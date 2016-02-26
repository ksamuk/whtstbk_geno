# per locus fst exploration

# libraries

install.packages("qqman")
library("qqman")
library("dplyr")
library("ggplot2")
library("tidyr")
library("ggthemes")

list.files("functions", full.names = TRUE) %>% sapply(.,source, verbose = FALSE, echo = FALSE) %>% invisible

# read in raw data
perloc_files <- list.files("data/stats", pattern = "50k", full.names = TRUE)

cbr_cmn_fst <- read.table(perloc_files[1], header = TRUE, stringsAsFactors = FALSE)
names(cbr_cmn_fst) <- c("chr", "pos1", "pos2", "n_snps", "w_fst", "fst")
cbr_cmn_fst$fst[cbr_cmn_fst$fst < 0] <- 0
cbr_cmn_fst <- cbr_cmn_fst[,c(1:3, 5)]

wht_cmn_fst <- read.table(perloc_files[3], header = TRUE, stringsAsFactors = FALSE)
names(wht_cmn_fst) <- c("chr", "pos1", "pos2", "n_snps", "w_fst", "fst")
wht_cmn_fst$fst[wht_cmn_fst$fst < 0] <- 0
wht_cmn_fst <- wht_cmn_fst[,c(1:3, 5)]

wht_cbr_fst <- read.table(perloc_files[2], header = TRUE, stringsAsFactors = FALSE)
names(wht_cbr_fst) <- c("chr", "pos1", "pos2", "n_snps", "w_fst", "fst")
wht_cbr_fst$fst[wht_cbr_fst$fst < 0] <- 0
wht_cbr_fst <- wht_cbr_fst[,c(1:3, 5)]

fst_df <- full_join(wht_cmn_fst, wht_cbr_fst, by = c("chr", "pos1", "pos2"))
fst_df <- full_join(fst_df, cbr_cmn_fst, by = c("chr", "pos1", "pos2"))
names(fst_df)[4:6] <- c("wht_cmn", "wht_cbr", "cbr_cmn")
fst_df$chr <- fst_df$chr %>% gsub("chr", "", .) %>% gsub("Un", "XXII", .) %>% as.roman %>% as.numeric


#### READ IN XTX FILES FROM BAYPASS

slugs <- c("wht_cmn", "wht_cbr", "cmn_cbr")

read_xtx_file <- function(slug){
  
  baypass_sites <- read.table("metadata/baypass_sites.txt", header = TRUE, stringsAsFactors = FALSE)
  
  # the baypass output files
  xtx_file <- list.files(paste0("data/baypass/", slug), pattern = paste0(slug, "_summary_pi_xtx"),full.names = TRUE)
  xtx_df <- read.table(xtx_file, header = TRUE, stringsAsFactors = FALSE)
  names(xtx_df) <- tolower(names(xtx_df)) 
  xtx_df <- cbind(baypass_sites, xtx_df[,-1])
  names(xtx_df)[7] <- paste0("xtx_", slug)
  xtx_df[,c(1:2, 7)]
}

xtx_list <- lapply(slugs, read_xtx_file)

xtx_df <- Reduce(function(x, y) merge(x, y, all=TRUE), xtx_list) %>%
  arrange(chr, pos)

#### Merge 

####

# join xtx files?
fx_df <- full_join(fst_df, xtx_df)

fx_df <- fst_df

wht_cmn_outlier <- is.outlier(fx_df$wht_cmn, cutoff = 0.99)
wht_cbr_outlier <- is.outlier(fx_df$wht_cbr, cutoff = 0.99)
cbr_cmn_outlier <- is.outlier(fx_df$cbr_cmn, cutoff = 0.99)

fx_df$wht_outlier <- (wht_cmn_outlier & wht_cbr_outlier) & !(cbr_cmn_outlier)

fst_long <- gather(fx_df, key = fst_type, value = fst, -chr, -pos1, -pos2, -wht_outlier)

fx_df %>%
  #sample_frac(0.01)%>%
  ggplot(aes(x = cbr_cmn, y = xtx_cmn_cbr)) +
  geom_point(size = 1)+
  facet_grid(fst_type~chr, scales = "free", space = "free_x", switch = "x") +
  theme_classic()+
  theme(axis.text.x = element_blank(), 
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        panel.margin = unit(0.1, "cm"),
        strip.background = element_blank())

fst_long %>%
  #filter(div_stat %in% c("wht_cmn", "wht_cbr", "cbr_cmn")) %>%
  #filter(chr %in% c(1,4,7)) %>%
  #sample_frac(0.01) %>%
  ggplot(aes(x = pos2, y = fst, color = wht_outlier)) +
  geom_point(size = 1)+
  facet_grid(fst_type~chr, scales = "free", space = "free_x", switch = "x") +
  theme_classic()+
  theme(axis.text.x = element_blank(), 
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        panel.margin = unit(0.1, "cm"),
        strip.background = element_blank())+
  scale_color_manual(values=c("black", "red"))

