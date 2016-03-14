# plot treemix output

library("dplyr")
list.files("functions", full.names = TRUE) %>% sapply(.,source, verbose = FALSE, echo = FALSE) %>% invisible


folder <- "data/treemix/wht_cmn_outgroups"

files <- list.files(folder, full.names = TRUE)
slugs <- strsplit(files, "\\.") %>% lapply(., function(x)x[1]) %>% unlist %>% unique

#plot_tree = function(stem, o = NA, cex = 1, disp = 0.003, plus = 0.01, flip = vector(), arrow = 0.05, scale = T, ybar = 0.1, mbar = T, plotmig = T, plotnames = T, xmin = 0, lwd = 1, font = 1)

# arrow = size of arrowhead
# lwd = lwd of tree (not arrows)
# font = ???
# disp = displacement of labels
# plus = x-stretch of tree
# o = read.table(o, as.is = T, comment.char = "", quote = "")

lapply(slugs[4], plot_tree, arrow = 0.1, lwd = 2, font = 1, disp = 0.003, plus = 0.01, arrow_lwd = 2, plotnames = FALSE, spoof_labels = c("a", "b"))

text(labels = c("a", "b", "c"), x = 0.005, y = -1, cex = 0.75)


pop_order <- data.frame(pop = c("AL_cmn","CL_cmn","CL_wht","DK_dk","GC_cbr","LN_cbr","MH_cmn","MH_wht","MR_cbr","RT_wht","SF_cmn","SF_wht", "SH_cmn","SH_wht","SK_cbr","SR_cmn","SR_wht"))
pop_order$cluster <- gsub("[^a-z]*", "", pop_order$pop)
pop_order <- pop_order %>%
  arrange(cluster)

lapply(slugs[2], plot_resid_fixed, pop_order = pop_order$pop)
