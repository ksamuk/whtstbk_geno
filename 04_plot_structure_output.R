
library(dplyr)
library(tidyr)
library(ggplot2)

list.files("functions", full.names = TRUE) %>% sapply(.,source, verbose = FALSE, echo = FALSE) %>% invisible

# plot structure output
#structure_files <- list.files("data/structure/results_rename", pattern = "simple.[0-9].meanQ",full.names = TRUE)
structure_files <- list.files("data/structure/results_full", pattern = "simple.[0-9].meanQ",full.names = TRUE)


# read in base metadata file
meta_df <- read.csv("metadata/sex_reg.csv")

# make df for each structure run
str_k2 <- read_structure(structure_files[2], meta_df)
str_k3 <- read_structure(structure_files[3], meta_df)
str_k4 <- read_structure(structure_files[4], meta_df)
str_k5 <- read_structure(structure_files[5], meta_df)
str_k6 <- read_structure(structure_files[6], meta_df)


str_k4 %>%
  arrange(pop)%>%
  #filter(!pop=="MM",!pop=="NHR",!pop=="QR",!pop=="MM", !pop=="LD", !pop=="CB", !pop=="SP") %>%
  ggplot(aes(x=id, y=q.value, fill=factor(k)))+
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
  facet_wrap(year~pop, scales = "free_x")

## Prepare meta data (as harvested from structure)

# species calls
sp_df <- str_k4 %>%
  filter(k == "k2") %>%
  mutate(species = ifelse(q.value > 0.01, "wht", "cmn"))%>%
  select(-q.value, -k)

# genetic sex calls
sex_df <- str_k4 %>%
  filter(k == "k3") %>%
  mutate(gen.sex = ifelse(q.value > 0.4, "F", "M")) %>%
  select(-q.value, -k)

# join species and sex
mega_meta <- left_join(sp_df, sex_df)

write.csv(mega_meta, "metadata/mega_meta.csv", quote = FALSE, row.names = FALSE)

