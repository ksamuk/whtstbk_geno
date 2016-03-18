# GO / QTL Anaylsis

library("dplyr")

# get gene list for all sites
sites_df <- read.table("metadata/all_whtstbk_sites.txt", h = T)
wind_df <- sites_df %>% 
  select(chr, window, fst_outlier_adj, xtx_outlier_adj) %>%
  distinct %>%
  arrange(chr, window)

# ENSEMBL gene annotations from glazer assembly
glazer_genes <- read.table("metadata/ensGene_revised.gtf", h = FALSE)

# subset and rename
glazer_genes <- glazer_genes[,c(1,3:5,10,12)]

# rename columns
names(glazer_genes) <- c("chr", "type", "pos1", "pos2", "gene_id", "trans_id")

glazer_genes <- filter(glazer_genes, chr != "chrM")

# chr to numeric
glazer_genes$chr <- glazer_genes$chr %>% gsub("chr", "", .) %>% gsub("Un", "XXII", .) %>% as.roman %>% as.numeric

# clean gene/transnames

glazer_genes$gene_id <- glazer_genes$gene_id %>% gsub(";", "", .)
glazer_genes$trans_id <- glazer_genes$trans_id %>% gsub(";", "", .)


glazer_genes <- glazer_genes %>%
  arrange(chr, pos1)
# join sites + genes

sites_df$pos

find_gene_by_chr_pos <- function(chr_target, pos_target, glazer_genes, slack = 0, window_size = 75000){
  
  if(window_size > 0){
    targ1 <- pos_target
    targ2 <- pos_target + (window_size - 2)
    
    gene_df1 <- glazer_genes %>%
      filter(chr == chr_target) %>%
      filter(pos1 - slack <= targ1, pos2 + slack >= targ2)
    
    gene_df2 <- glazer_genes %>%
      filter(chr == chr_target) %>%
      filter(pos1 - slack <= targ2, pos2 + slack >= targ2)
    
    unique(gene_df1$gene_id, gene_df2$gene_id)
    
  } else{
    
    gene_df <- glazer_genes %>%
      filter(chr == chr_target) %>%
      filter(pos1 - slack <= pos_target, pos2 + slack >= pos_target)
    
    unique(gene_df$gene_id)
    
  }
  
}

tmp_wind <- sites_df %>% filter(chr == 1) %>% select(window) %>% unlist

find_gene_by_chr_pos(1, tmp_wind[50], glazer_genes)
