# plot treemix output

library("dplyr")
list.files("functions", full.names = TRUE) %>% sapply(.,source, verbose = FALSE, echo = FALSE) %>% invisible


folder <- "data/treemix/output/wht_cmn_outgroups"

files <- list.files(folder, full.names = TRUE)
slugs <- strsplit(files, "\\.") %>% lapply(., function(x)x[1]) %>% unlist %>% unique

lapply(slugs[13], plot_tree)
#text(labels = slugs[1], x = 0.005, y = -1, cex = 0.75)


pop_order <- data.frame(pop = c("AL_cmn","CL_cmn","CL_wht","DK_dk","GC_cbr","LN_cbr","MH_cmn","MH_wht","MR_cbr","RT_wht","SF_cmn","SF_wht", "SH_cmn","SH_wht","SK_cbr","SR_cmn","SR_wht"))
pop_order$cluster <- gsub("[^a-z]*", "", pop_order$pop)
pop_order <- pop_order %>%
  arrange(cluster)

lapply(slugs[3], plot_resid_fixed, pop_order = pop_order$pop)
