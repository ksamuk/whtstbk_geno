library("ggplot2")
library("dplyr")

fst_dist <- read.table("metadata/fst_distances.txt", h = T)

fst_dist %>%
  #filter(type %in% c("wht_wht", "cmn_cmn", "wht_cmn", :))%>%
  mutate(coast_distance = least.cost.distance + dist.to.coast1 + dist.to.coast2) %>%
  #mutate(comparison = ifelse(type %in% c("wht_wht", "cmn_cmn", "cbr_cbr"), "same", "different")) %>%
  ggplot(aes(x = log(least.cost.distance), y = fst, color = geography)) +
  geom_point()+
  #facet_wrap(~type) +
  geom_smooth(method = "lm", se = FALSE)

fst_dist %>%
  #filter(type %in% c("wht_wht", "cmn_cmn", "wht_cmn", :))%>%
  mutate(coast_distance = least.cost.distance + dist.to.coast1 + dist.to.coast2) %>%
  #mutate(comparison = ifelse(type %in% c("wht_wht", "cmn_cmn", "cbr_cbr"), "same", "different")) %>%
  ggplot(aes(x = geography, y = fst, color = type)) +
  geom_boxplot()
