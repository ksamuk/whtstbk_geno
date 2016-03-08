
library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)
library(wesanderson)
library(grid)

list.files("functions", full.names = TRUE) %>% sapply(.,source, verbose = FALSE, echo = FALSE) %>% invisible

# plot structure output
#structure_files <- list.files("data/structure/results_rename", pattern = "simple.[0-9].meanQ",full.names = TRUE)
structure_files <- list.files("data/structure/results_2014", pattern = "simple.*meanQ",full.names = TRUE)

# read in base metadata file
meta_df <- read.csv("metadata/mega_meta.csv")
structure_ids <- read.table("data/structure/structure_ids.txt")

# make df for each structure run
str_k2 <- read_structure(structure_files[2], structure_ids, meta_df)
str_k3 <- read_structure(structure_files[3], structure_ids, meta_df)
str_k4 <- read_structure(structure_files[4], structure_ids, meta_df)
str_k5 <- read_structure(structure_files[5], structure_ids, meta_df)
#str_k6 <- read_structure(structure_files[6], meta_df)

# DISTRUCT-esque plots

# if facetting, this looks pray good bruh
str_k3 %>%
  filter(!is.na(cluster)) %>%
  group_by(id) %>%
  mutate(major_k_qval = max(q.value)) %>%
  mutate(major_k = k[q.value == max(q.value)]) %>% 
  ungroup %>%
  arrange(major_k, major_k_qval) %>%
  ungroup %>%
  mutate(id = factor(id, levels = as.character(id))) %>%
  ggplot(aes(x = id, y = q.value, fill = factor(k)))+
  geom_bar(stat = "identity", width = 1)+
  facet_grid(~region, scales = "free", space = "free", switch = "both") +
  theme_classic()+
  theme(axis.text = element_blank(), 
        axis.ticks = element_blank(), 
        axis.line = element_blank(),
        axis.title = element_blank(),
        strip.text = element_text(size=12),
        panel.margin = unit(0.1, "lines"), 
        legend.position = "none", 
        strip.background = element_blank())


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

