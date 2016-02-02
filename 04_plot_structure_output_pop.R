
library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)
library(wesanderson)
library(grid)

list.files("functions", full.names = TRUE) %>% sapply(.,source, verbose = FALSE, echo = FALSE) %>% invisible

# plot structure output
#structure_files <- list.files("data/structure/results_rename", pattern = "simple.[0-9].meanQ",full.names = TRUE)
structure_files <- list.files("data/structure/results_pop", pattern = "meanQ",full.names = TRUE)

# read in base metadata file
meta_df <- read.csv("metadata/mega_meta.csv")

tmp <- read_structure_pop(structure_files[12], meta_df)

k1_vals <- str_k2 %>%
  filter(year != 2012) %>%
  filter(k == "k1") %>%
  select(id, q.value) %>%
  rename(k1_val = q.value )

k3_vals <- str_k4 %>%
  filter(year != 2012) %>%
  filter(k == "k3") %>%
  select(id, q.value) %>%
  rename(k3_val = q.value )

tmp <- left_join(str_k2, k1_vals)
tmp <- left_join(tmp, k3_vals)


tmp %>%
  ggplot(aes(x = id, y = q.value, fill = factor(k)))+
  geom_bar(stat = "identity", width = 1, color = "black")+
  theme_classic()+
  facet_wrap(pop~gen.sex)+
	scale_fill_brewer(palette = "Set1")

tmp %>%
  filter(year != 2012) %>%
  group_by(id) %>%
  mutate(max_q = max(q.value)) %>% 
  mutate(major_k = k[which(q.value == max_q)]) %>%
  ungroup()%>%
  mutate(pop = ifelse(grepl("CP", pop), "CP", as.character(pop))) %>%
  mutate(id = reorder(id, as.numeric(k1_val))) %>%
  filter(!pop=="MM",!pop=="NHR",!pop=="QR",!pop=="MM", !pop=="LD", !pop=="CB", !pop=="SP") %>%
  #mutate(id = reorder(id, as.numeric(major_k))) %>%
  ungroup() %>%
  mutate(id = reorder(id, as.numeric(k3_val))) %>%
  #filter(!pop=="MM",!pop=="NHR",!pop=="QR",!pop=="MM", !pop=="LD", !pop=="CB", !pop=="SP") %>%
  ggplot(aes(x = id, y = q.value, fill = factor(k)))+
  #ggplot(aes(x=id, y=q.value, fill=factor(pop)))+
  geom_bar(stat="identity", width = 2, color = "black")+
  #geom_bar(aes(x=id,fill=pop, y=.01),stat="identity",width=1,position="stack")+
  #geom_text(aes(label=pop))+
  theme_classic()+
  theme(axis.text = element_blank(), 
        axis.ticks = element_blank(), 
        axis.line = element_blank(),
        axis.title = element_blank(),
        strip.text = element_text(size=12),
        panel.margin = unit(0.1, "lines"))+
  facet_wrap(pop~gen.sex, scales = "free_x")+
  scale_fill_brewer(palette = "Set1")

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

