
library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)
library(wesanderson)

list.files("functions", full.names = TRUE) %>% sapply(.,source, verbose = FALSE, echo = FALSE) %>% invisible

# plot structure output
#structure_files <- list.files("data/structure/results_rename", pattern = "simple.[0-9].meanQ",full.names = TRUE)
structure_files <- list.files("data/structure/results_full", pattern = "simple.[0-9].meanQ",full.names = TRUE)


# read in base metadata file
meta_df <- read.csv("metadata/mega_meta.csv")

# make df for each structure run
str_k2 <- read_structure(structure_files[2], meta_df)
str_k3 <- read_structure(structure_files[3], meta_df)
str_k4 <- read_structure(structure_files[4], meta_df)
str_k5 <- read_structure(structure_files[5], meta_df)
str_k6 <- read_structure(structure_files[6], meta_df)

k3_vals <- str_k4 %>%
  filter(year != 2012) %>%
  filter(k == "k3") %>%
  select(id, q.value) %>%
  rename(k3_val = q.value )

tmp <- left_join(str_k4, k3_vals)


tmp %>%
	filter(year != 2012) %>%
  group_by(sex) %>%
  mutate(max_q = max(q.value)) %>% 
  mutate(major_k = k[which(q.value == max_q)]) %>%
  #mutate(id = reorder(id, as.numeric(major_k))) %>%
  mutate(id = reorder(id, as.numeric(k3_val))) %>%
  #filter(!pop=="MM",!pop=="NHR",!pop=="QR",!pop=="MM", !pop=="LD", !pop=="CB", !pop=="SP") %>%
  ggplot(aes(x = id, y = q.value, fill = factor(k)))+
  #ggplot(aes(x=id, y=q.value, fill=factor(pop)))+
  geom_bar(stat="identity", width=1, color = "black")+
  #geom_bar(aes(x=id,fill=pop, y=.01),stat="identity",width=1,position="stack")+
  #geom_text(aes(label=pop))+
  theme_classic()+
  theme(axis.text=element_blank(), 
        axis.ticks=element_blank(), 
        axis.line=element_blank(),
        axis.title=element_blank(),
        strip.text=element_text(size=12))+
  facet_wrap(~gen.sex, scales = "free_x")+
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

