################################################################################
# initials
################################################################################

library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)
library(wesanderson)
library(grid)
library(cowplot)
library(ggthemes)

list.files("functions", full.names = TRUE) %>% sapply(.,source, verbose = FALSE, echo = FALSE) %>% invisible

################################################################################
# raw data
################################################################################

# plot structure output
#structure_files <- list.files("data/structure/results_rename", pattern = "simple.[0-9].meanQ",full.names = TRUE)
structure_files <- list.files("data/structure/results_2014", pattern = "simple.*meanQ",full.names = TRUE)

# read in base metadata file
meta_df <- read.csv("metadata/mega_meta.csv")

# long region names
region_short <- c("AN", "CB", "HA", "GY")
region_names <- c("Antigonish", "Bras d'Or", "Halifax", "Guysborough")

# meta data
meta_df <- meta_df %>%
mutate(region = region_names[match(region, region_short)])

# separated structure ids
structure_ids <- read.table("data/structure/structure_ids.txt")

################################################################################
# read in structure data (k1-k3 are significant)
################################################################################

# make df for each structure run
str_k2 <- read_structure(structure_files[2], structure_ids, meta_df)
str_k3 <- read_structure(structure_files[3], structure_ids, meta_df)

# color palatte
palette4 <- clpalette('694737')$colors %>% as.character %>% paste0("#", .) %>% rev

# DISTRUCT-esque plots

################################################################################
# my personal take on a distruct plot, using ggplot
################################################################################

myPalette <- colorRampPalette(rev(brewer.pal(3, "Spectral")), space="Lab")

plot_structure <- function(x, k, labels = TRUE){
  
  if(labels){
    label_theme <- theme(strip.text = element_text(size=16),
                         plot.margin = unit(c(0,0.5,0.5,0.25), "cm"))
  } else{
    label_theme <- theme(strip.text = element_blank(),
                         plot.margin = unit(c(0,0.5,0.0,0.25), "cm"))
  }
  
  k_val <- paste0("K = ", k)
  #y_lab <- paste0("K = ", k, "\n", "q-value")
  #y_lab <- paste0("q-value, ", "K = ", k)
  
  ann_text <- data.frame(region = "Antigonish", lab = k_val, q.value = 0.8, id = "NU1", k = "k1")
  
  x %>%
  filter(!is.na(cluster)) %>%
    group_by(id) %>%
    mutate(major_k_qval = max(q.value)) %>%
    mutate(major_k = k[q.value == max(q.value)]) %>% 
    ungroup %>%
    arrange(major_k, major_k_qval) %>%
    ungroup %>%
    mutate(id = factor(id, levels = as.character(id))) %>%
    ggplot(aes(x = id, y = q.value, fill = factor(k)))+
    geom_bar(stat = "identity", width = 1.1)+
    facet_grid(~region, scales = "free", space = "free", switch = "x") +
    theme_classic()+
    theme(axis.text.x = element_blank(), 
          axis.ticks = element_blank(), 
          axis.line = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_text(size=16),
          panel.margin = unit(0.1, "lines"), 
          legend.position = "none", 
          strip.background = element_blank(),
          strip.switch.pad.grid = unit(0.0, "lines"))+
    label_theme+
    geom_text(data = ann_text, label = k_val, aes(y = 0.96, x = 8, fontface = "bold"), size = 5)+
    ylab("q-value")+
    scale_fill_discrete()
}

################################################################################
# plot k2-3
################################################################################

# if facetting, this looks pray good bruh
k2 <- str_k2 %>%
  plot_structure(labels = FALSE, k = 2)

k3 <- str_k3 %>%
  plot_structure(k = 3)

figure_3 <- plot_grid(k2, k3, ncol = 1, rel_heights = c(0.85,1))
ggsave("figures/Figure3.pdf", plot = figure_3, width = 11, height = 8.5)

################################################################################
# not run 
################################################################################

# if k = 2 and NOT facetting, this looks better
str_k2 %>%
  filter(!is.na(cluster)) %>%
  group_by(id) %>%
  mutate(k1val = q.value[k=="k1"]) %>%
  ungroup %>%
  mutate(id = reorder(id, -k1val)) %>%
  ggplot(aes(x = id, y = q.value, fill = factor(k)))+
  geom_bar(stat = "identity", width = 1)+
  #facet_grid(~major_k, scales = "free_x", switch = "x") +
  theme_classic()+
  theme(axis.text = element_blank(), 
        axis.ticks = element_blank(), 
        axis.line = element_blank(),
        axis.title = element_blank(),
        strip.text = element_text(size=12),
        panel.margin = unit(0.1, "lines"), 
        legend.position = "none", 
        strip.background = element_blank())

## Prepare meta data (as harvested from structure)

# species calls
sp_df <- str_k4 %>%
  filter(k == "k2") %>%
  mutate(species = ifelse(q.value > 0.01, "wht", "cmn"))%>%
  select(-q.value, -k)

# genetic sex calls
k4_spread <- str_k4 %>% spread(key = k, value = q.value)

sex_df <- k4_spread %>%
  mutate(gen.sex = ifelse(k3 + k2 > 0.5, "F", "M")) %>%
  select(-k1, -k2, -k3, -k4)

# join species and sex
mega_meta <- left_join(sp_df, sex_df)

#write.csv(mega_meta, "metadata/mega_meta.csv", quote = FALSE, row.names = FALSE)

