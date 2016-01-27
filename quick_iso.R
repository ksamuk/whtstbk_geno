# quick isotope analysis

library(ggplot2)
library(dplyr)

iso <- read.csv("isotope_data.csv")
struc <- read.table("whtcmn_faststructure.txt", header = TRUE, stringsAsFactors = FALSE)

struc <- struc %>%
	filter(k.value.run == 2) %>%
	select(id, membership)

iso <- iso %>%
	mutate(pop = gsub("[^A-Z]*", "", as.character(id))) %>%
	mutate(cn.ratio = C.amount / N.amount)

iso <- left_join(iso, struc)

#iso$d15N.resid <- lm(data = iso, d15N ~ pop + cn.ratio) %>% residuals
#iso$d13C.resid <- lm(data = iso, d13C ~ pop + cn.ratio) %>% residuals

iso$d15N.resid <- lm(data = iso, d15N ~ pop + cn.ratio) %>% residuals
iso$d13C.resid <- lm(data = iso, d13C ~ pop + cn.ratio ) %>% residuals

iso %>%
	ggplot(aes(x = d13C.resid, y = d15N.resid, label = id, color = membership)) +
	geom_point(size = 3) +
	facet_wrap(~pop)
