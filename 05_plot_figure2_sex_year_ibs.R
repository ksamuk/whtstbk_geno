################################################################################
# initials
################################################################################

# libraries
library("ggplot2")
library("gdsfmt")
library("SNPRelate")
library("dplyr")
library("tidyr")
library("ggthemes")
library("RColorBrewer")
library("bigmemory")
library("cowplot")

list.files("functions", full.names = TRUE) %>% sapply(.,source, verbose = FALSE, echo = FALSE) %>% invisible

select <- dplyr::select

################################################################################
# initials
################################################################################

# pca and ibs data files
pca_df <- read.table("metadata/pca_df.txt", h = T)
diss_df <- read.table("metadata/diss_df.txt", h = T)
var_prop <- read.table("metadata/pca_var_df.txt", h = T)[,1]

################################################################################
# plots
################################################################################
# lost all the code here by accident...sad moment, didn't commit. 

year_pca <- pca_df %>%
  ggplot(aes(x = EV1, y = EV2, fill = as.factor(year)))+
  geom_point(pch = 21, size = 4, color = "black")+
  xlab(paste0("PC1 ", "(",var_prop[1], "%)"))+
  ylab(paste0("PC2 ", "(",var_prop[2], "%)"))+
  theme_base() +
  theme(legend.background = element_blank(),
        text = element_text(size = 14),
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        legend.justification = c(1, 1), 
        legend.position = c(1, 1),
        legend.position = c(1, 1),
        #axis.title.y = element_text(margin=margin(0,6,0,0), face = "bold"),
        #axis.title.x = element_text(margin=margin(12,0,0,0), face = "bold"),
        plot.background = element_rect(color = NA))+
  scale_fill_fivethirtyeight()

sex_pca <- pca_df %>%
  mutate(sex = ifelse(sex == "M", "Male", "Female")) %>%
  filter(year == 2014) %>%
  ggplot(aes(x = EV1, y = EV2, fill = as.factor(sex)))+
  geom_point(pch = 21, size = 4, color = "black")+
  xlab(paste0("PC1 ", "(",var_prop[1], "%)"))+
  ylab(paste0("PC2 ", "(",var_prop[2], "%)"))+
  theme_base() +
  theme(legend.background = element_blank(),
        text = element_text(size = 14),
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        legend.justification = c(1, 1), 
        legend.position = c(1, 1),
        legend.position = c(1, 1),
        #axis.title.y = element_text(margin=margin(0,6,0,0), face = "bold"),
        #axis.title.x = element_text(margin=margin(12,0,0,0), face = "bold"),
        plot.background = element_rect(color = NA))+
  scale_fill_fivethirtyeight()

year_ibs <- diss_df %>%
  #filter(diss < 0.79) %>%
  group_by(id, year_type, cluster_type) %>%
  summarise(mean_diss = mean(diss)) %>% ungroup %>%
  mutate(year_type = ifelse(year_type == "different", "Different\nYear", "Same\nYear")) %>%
  mutate(cluster_type = ifelse(cluster_type == "different", "Between Cluster", "Within Cluster")) %>%
  ggplot(aes(x = year_type, y = mean_diss, fill = as.factor(cluster_type)))+
  #geom_boxplot()+
  geom_jitter(pch = 21, size = 2, color = "black")+
  stat_summary(fun.y=median, fun.ymin=median, fun.ymax=median, 
               geom="crossbar", width=0.9, color = "black", fatten = 6)+
  facet_grid(.~cluster_type, switch = "x")+
  ylab("% Shared SNPs")+
  theme_base() +
  theme(legend.background = element_blank(),
        text = element_text(size = 14),
        axis.title.x = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        legend.justification = c(1, 1), 
        legend.position = "none",
        #axis.title.y = element_text(margin=margin(0,6,0,0), face = "bold"),
        #axis.title.x = element_text(margin=margin(12,0,0,0), face = "bold"),
        plot.background = element_rect(color = NA))+
  scale_fill_fivethirtyeight()

sex_ibs <- diss_df %>%
  #filter(diss < 0.79) %>%
  filter(year_1 == 2014 & year_2 ==2014) %>%
  filter(sex_1 == sex_2) %>%
  mutate(sex_type = ifelse(sex_1 == "F", "Females", "Males")) %>%
  #mutate(cluster_type = ifelse(cluster_1 == "wht" & cluster_2 == "wht", "White vs.\n White", ifelse(cluster_1 == "cmn" & cluster_2 == "cmn", "Common vs.\n Common", "White vs.\n Common"))) %>%
  mutate(cluster_type = ifelse(cluster_type == "different", "Between Cluster", "Within Cluster")) %>%
  group_by(cluster_type,sex_type, id) %>%
  summarise(mean_diss = mean(diss)) %>% ungroup %>%
  ggplot(aes(x = sex_type, y = mean_diss, fill = as.factor(cluster_type)))+
  #geom_boxplot()+
  geom_jitter(pch = 21, size = 2, color = "black")+
  stat_summary(fun.y=median, fun.ymin=median, fun.ymax=median, 
               geom="crossbar", width=0.9, color = "black", fatten = 6)+
  facet_grid(.~cluster_type, switch = "x")+
  ylab("% Shared SNPs")+
  facet_grid(.~cluster_type, switch = "x")+
  theme_base() +
  theme(legend.background = element_blank(),
        text = element_text(size = 14),
        axis.title.x = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        legend.justification = c(1, 1), 
        legend.position = "none",
        #axis.title.y = element_text(margin=margin(0,6,0,0), face = "bold"),
        #axis.title.x = element_text(margin=margin(12,0,0,0), face = "bold"),
        plot.background = element_rect(color = NA))+
  scale_fill_fivethirtyeight()
  
  
figure_2 <- plot_grid(year_pca,year_ibs, sex_pca, sex_ibs, scale = 0.95)  

ggsave(figure_2, file = "figures/figure2.pdf", height = 8.5, width = 8.5)
 

################################################################################
# LN probe
################################################################################

diss_df %>%
  filter(cluster_1 != "cmn" & cluster_2 != "cmn") %>%
  filter(grepl("CL", id)&grepl("CL", id2)) %>% View
mutate(is_ln_20 = (id == "CL50" | id2 == "CL50")) %>%
  ggplot(aes(x = cluster_type, y = diss, color = is_ln_20))+
  geom_boxplot()