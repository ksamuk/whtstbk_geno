# plot allele freq data


library("stringr")
library("dplyr")
library("tidyr")
library("readr")
library(ggplot2)
library(hexbin)
library(viridis)

list.files("functions", full.names = TRUE) %>% sapply(.,source, verbose = FALSE, echo = FALSE) %>% invisible

wht_frq <- read_delim("data/allele_freq_by_pop.txt", delim = " ")
wht_fst <- read_delim("data/fst_outlier_list.txt", delim = " ")

wht_fst %>%
  sample_frac(0.5) %>%
  #filter(fst.outlier.male == TRUE | fst.outlier.female == TRUE) %>%
  ggplot(aes(x = pos, y = as.numeric(fst.male), color = fst.outlier.male)) +
  geom_point() +
  geom_point(aes(x = pos, y = as.numeric(fst.female), color = fst.outlier.female)) +
  facet_wrap(~chr)

chr4_frq <- wht_frq %>%
  filter(chr == 4) %>%
  #filter(pos > 2000000) %>%
  arrange(pop, chr, pos, allele) %>%
  filter(!is.na(allele))
  
chr4_frq <- chr4_frq %>%
  group_by(pop, chr, pos) %>%
  summarise(alt_count = sum(n[allele == 2]*2, n[allele == 1], na.rm = TRUE), total_count = sum(n, na.rm = TRUE)*2) %>%
  mutate(alt_frq = alt_count / total_count) %>%
  mutate(species = ifelse(grepl("wht", pop), "white", "common")) %>%
  mutate(year = ifelse(grepl("2012", pop), 2012, 2014))

chr4_frq %>%
  ungroup %>%
  filter(year > 2012) %>%
  group_by(chr, pos, species, pop) %>%
  summarise(alt_frq = mean(alt_frq)) %>%
  ungroup %>%
  #mutate(pos = as.numeric(as.factor(pos))) %>%
  filter(!is.na(alt_frq)) %>%
  filter(!grepl("RT_cmn|GC_wht", pop)) %>%
  mutate(pop = reorder(pop, as.numeric(as.factor(species)))) %>%
  #mutate(species = as.numeric(as.factor(species))) %>%
  #mutate(pop = reorder(pop, species)) %>%
  ggplot(aes(x = pos, y = pop, color = alt_frq))+
  geom_tile()+
  scale_color_viridis()
  
