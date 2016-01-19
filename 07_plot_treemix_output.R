# plot treemix output

library("dplyr")
list.files("functions", full.names = TRUE) %>% sapply(.,source, verbose = FALSE, echo = FALSE) %>% invisible


folder <- "data/treemix/output/noAL_rootCP"

files <- list.files(folder, full.names = TRUE)
slugs <- strsplit(files, "\\.") %>% lapply(., function(x)x[1]) %>% unlist %>% unique

lapply(slugs, plot_tree)

lapply(slugs, plot_resid)
