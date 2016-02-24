# per locus fst exploration

# libraries

install.packages("qqman")
library("qqman")
library("dplyr")
library("ggplot2")
library("tidyr")
library("ggthemes")
# read in raw data
perloc_files <- list.files("data/stats", pattern = "perloc", full.names = TRUE)

cbr_cmn_fst <- read.table(perloc_files[1], header = TRUE, stringsAsFactors = FALSE)
names(cbr_cmn_fst) <- c("chr", "pos", "fst")
cbr_cmn_fst$fst[cbr_cmn_fst$fst < 0] <- 0

wht_cmn_fst <- read.table(perloc_files[2], header = TRUE, stringsAsFactors = FALSE)
names(wht_cmn_fst) <- c("chr", "pos", "fst")
wht_cmn_fst$fst[wht_cmn_fst$fst < 0] <- 0


wht_cbr_fst <- read.table(perloc_files[3], header = TRUE, stringsAsFactors = FALSE)
names(wht_cbr_fst) <- c("chr", "pos", "fst")
wht_cbr_fst$fst[wht_cbr_fst$fst < 0] <- 0

fst_df <- full_join(wht_cmn_fst, wht_cbr_fst, by = c("chr", "pos"))
fst_df <- full_join(fst_df, cbr_cmn_fst, by = c("chr", "pos"))
names(fst_df)[3:5] <- c("wht_cmn", "wht_cbr", "cbr_cmn")
fst_df$chr <- fst_df$chr %>% gsub("chr", "", .) %>% gsub("Un", "XXII", .) %>% as.roman %>% as.numeric


# join xtx files?
fx_df <- left_join(fst_df, xtx_df)

fst_long <- gather(fst_df, key = fst_type, value = fst, -chr, -pos)

fx_df %>%
  #sample_frac(0.01)%>%
  ggplot(aes(x = m_xtx, y = wht_cbr)) +
  geom_point(size = 1)+
  facet_grid(fst_type~chr, scales = "free", space = "free_x", switch = "x") +
  theme_classic()+
  theme(axis.text.x = element_blank(), 
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        panel.margin = unit(0.1, "cm"),
        strip.background = element_blank())