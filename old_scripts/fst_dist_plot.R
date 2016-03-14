library("ggplot2")
library("dplyr")

fst_dist <- read.table("metadata/fst_distances.txt", h = T)

fst_dist %>%
  filter(type %in% c("wht_wht", "cmn_cmn", "wht_cmn"))%>%
  mutate(coast_distance = least.cost.distance + dist.to.coast1 + dist.to.coast2) %>%
  mutate(comparison = ifelse(type %in% c("wht_wht", "cmn_cmn", "cbr_cbr"), "same", "different")) %>%
  ggplot(aes(x = coast_distance, y = fst, color = type)) +
  geom_point()+
  geom_smooth(method = "lm")