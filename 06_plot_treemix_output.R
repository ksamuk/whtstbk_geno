 # plot treemix output

################################################################################
# Libraries
################################################################################

library("dplyr")
library("wesanderson")
library("ape")
library("phangorn")
list.files("functions", full.names = TRUE) %>% sapply(.,source, verbose = FALSE, echo = FALSE) %>% invisible

# treemix trees with little campbell river (pacific) + denmark outgroups
lcr_folder <- "data/treemix/lcr_dk_outgroups"
lcr_files <- list.files(lcr_folder, full.names = TRUE)
lcr_slugs <- strsplit(lcr_files, "\\.") %>% lapply(., function(x)x[1]) %>% unlist %>% unique

# treemix trees with denmark (atlantic) outgroups only
dk_folder <- "data/treemix/wht_cmn_outgroups"
dk_files <- list.files(dk_folder, full.names = TRUE)
dk_slugs <- strsplit(dk_files, "\\.") %>% lapply(., function(x)x[1]) %>% unlist %>% unique

################################################################################
# set up consensus tree
################################################################################

# bootstrap trees with denmark (atlantic) outgroups only
boot_file <- "data/treemix/whtstbk_2014cat_trees.tre"
boot_tree <- read.tree(boot_file)

tree_tmp <- consensus(boot_tree, p = 0.5)
prop.clades(tree_tmp, boot_tree)->tree_tmp$node.label


################################################################################
# plot treemix wrapper
################################################################################

#plot_tree = function(stem, o = NA, cex = 1, disp = 0.003, plus = 0.01, flip = vector(), arrow = 0.05, scale = T, ybar = 0.1, mbar = T, plotmig = T, plotnames = T, xmin = 0, lwd = 1, font = 1)

# arrow = size of arrowhead
# lwd = lwd of tree edges (not arrows)
# NEW arrow_lwd = lwd of migration edges
# font = ???
# disp = displacement of labels
# plus = x-stretch of tree

# hack to access labels and colors directly

plot_treemix_wrapper <- function(slug, ...){
  
suppressWarnings(spoof_labels <- get_treemix_d_object(slug))
spoof_labels <- spoof_labels[,2][!is.na(spoof_labels[,2])]
spoof_cols <- spoof_labels %>% gsub("[A-Z_]*", "", .) %>% as.factor %>% as.numeric

spoof_cols <- rev(brewer.pal(length(unique(spoof_cols)), "Set1"))[spoof_cols]

spoof_labels <- spoof_labels %>% gsub("[a-z_]*", "", .) %>% gsub("NA", "", .) 

suppressWarnings(plot_tree(slug, spoof_labels = spoof_labels, spoof_cols = spoof_cols, ...))
}

################################################################################
# plot figure 4
################################################################################

# base lcr/dk tree
plot_treemix_wrapper(lcr_slugs[23], arrow = 0.1, lwd = 3, font = 2, disp = 0.0005, 
                     plus = 0.005, arrow_lwd = 3, plotnames = TRUE, cex = 1, use_viridis = FALSE, use_alpha = FALSE,
                     arrow_lty = 1, ybar = 0.5, shadow = 0.1)

pdf(file = "figures/Figure4.pdf", width = 8.5, height = 4.25)

margins_trees <- c(4,0,0,0)
par(mfrow=c(1,2), mar = margins_trees)

# consensus tree
plot(tree_tmp, no.margin = FALSE)
nodelabels(tree_tmp$node.label/1000, frame = "none", adj = c(0.5, 0))

# best mig edge tree
plot_treemix_wrapper(dk_slugs[4], arrow = 0.2, lwd = 3, font = 2, disp = 0.0015, 
                     plus = 0.005, arrow_lwd = 3, plotnames = TRUE, cex = 1, use_viridis = FALSE, use_alpha = FALSE,
                     arrow_lty = 1, ybar = 0.5, shadow = 0.05)
dev.off()

################################################################################
# plot figure s2
################################################################################


pdf(file = "figures/FigureS2.pdf", width = 4.25, height = 4.25)

plot_treemix_wrapper(lcr_slugs[23], arrow = 0.2, lwd = 3, font = 2, disp = 0.0015, 
                     plus = 0.005, arrow_lwd = 3, plotnames = TRUE, cex = 1, use_viridis = FALSE, use_alpha = FALSE,
                     arrow_lty = 1, ybar = -10, shadow = 0.05)

dev.off()

# get mean branch lengths for nova scotia to pacific split

boot_file <- "data/treemix/outgroup_bootstrapcat_trees.tre"
boot_tree <- read.tree(boot_file)


compute_branch_length <- function(tree, node1 = 23, node2 = 21){
  
  dist.nodes(tree)[node1, node2]
}

atl_pac_branch <- lapply(boot_tree, compute_branch_length) %>% unlist

atl_mean <- lm(atl_pac_branch~1)
# mean: 0.04992031

confint(atl_mean)
#2.5 %     97.5 %
#  (Intercept) 0.04986507 0.04997555

compute_wht_cmn_length <- function(tree, node1 = 23, node2 = 21){
  
  dist.nodes(tree)[node1, node2]
}

# ortis et al:




################################################################################
# plot residuals
################################################################################

pop_order <- data.frame(pop = c("DK_dk", "LC_NA", "AL_cmn","CL_cmn","CL_wht","GC_cbr","LN_cbr","MH_cmn","MH_wht","MR_cbr","RT_wht","SF_cmn","SF_wht", "SH_cmn","SH_wht","SK_cbr","SR_cmn","SR_wht", "CP_cmn"))
pop_order$cluster <- gsub("[^a-z]*", "", pop_order$pop)
pop_order <- pop_order %>%
  arrange(cluster)

lapply(slugs[11], plot_resid_fixed, pop_order = as.character(pop_order$pop))

plot_resid(slugs[12])

################################################################################
# harvest migration edge pvalues
################################################################################

# this ain't gonna be pretty

treeout_files <- grep("treeout", dk_files, value = TRUE) %>% grep("m0", ., invert = TRUE, value = TRUE)

extract_mig_pvals <- function(tree_out_file){
  
  num_mig <- tree_out_file %>% gregexpr("m[0-9]", .) %>% as.numeric
  num_mig <- substr(tree_out_file, num_mig, num_mig+1) %>% gsub("m", "",.)
  
  treeout_raw <- scan(tree_out_file, what = "character")
  
  treeout_raw <- treeout_raw[-1] 
  
  n_mig <- (treeout_raw[-1] %>% length)/6
  
  indexes <- split(1:length(treeout_raw), ceiling(seq_along(1:length(treeout_raw))/6))
  
  mig_df <- lapply(indexes, function(x) data.frame(t(treeout_raw[x])) )
  mig_df <- bind_rows(mig_df)
  names(mig_df) <- c("mig_w", "mig_w_jack", "mig_se", "pval", "source", "sink")
  cbind(num_mig, mig_df)
}

mig_df <- lapply(treeout_files[1:6], extract_mig_pvals)
mig_df <- bind_rows(mig_df)

################################################################################
# harvest likelihoods
################################################################################

llik_files <- grep("llik", dk_files, value = TRUE)

extract_llik <- function(llik_file){
  
  num_mig <- llik_file %>% gregexpr("m[0-9]", .) %>% as.numeric
  num_mig <- substr(llik_file, num_mig+1, num_mig+1) %>% as.numeric
  
  llik_raw <- scan(llik_file, what = "character", sep ="\n")
  
  llik <- llik_raw[2] %>% strsplit(split = ":") %>% lapply(function(x)x[2]) %>% gsub(" ", "", .) %>% as.numeric
  
  data.frame(num_mig, llik)
}

llik_df <- lapply(llik_files[1:6], extract_llik)
llik_df <- bind_rows(llik_df)

lik_ratio_pval <- list()
lik_ratio_pval[1] <- NA

lik_ratio_chisq <- list()
lik_ratio_chisq[1] <- NA

for (i in 2:nrow(llik_df)){
  chi_sq <- (llik_df$llik[i] - llik_df$llik[i-1])*2
  lik_ratio_chisq[i] <- chi_sq
  lik_ratio_pval[i] <- pchisq(chi_sq, df = 1, lower.tail = FALSE)
  
}

data.frame(llik_df, unlist(lik_ratio_chisq), unlist(lik_ratio_pval))


