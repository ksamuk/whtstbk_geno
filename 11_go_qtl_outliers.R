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

# no mitochondria
glazer_genes <- filter(glazer_genes, chr != "chrM")

# chr to numeric
glazer_genes$chr <- glazer_genes$chr %>% gsub("chr", "", .) %>% gsub("Un", "XXII", .) %>% as.roman %>% as.numeric

# clean gene/transnames
glazer_genes$gene_id <- glazer_genes$gene_id %>% gsub(";", "", .)
glazer_genes$trans_id <- glazer_genes$trans_id %>% gsub(";", "", .)

# arrange
glazer_genes <- glazer_genes %>%
  arrange(chr, pos1)

# join sites + genes
sites_df$pos

find_gene_by_chr_pos <- function(index, glazer_genes, slack = 0, window_size = 75000){
  
  if(window_size > 0){
    
    chr_target <- sites_df$chr[index]
    pos_target <- sites_df$window[index]
    
    targ1 <- pos_target
    targ2 <- pos_target + (window_size - 2)
    
    gene_df1 <- glazer_genes %>%
      filter(chr == chr_target) %>%
      filter(pos1 - slack <= targ1, pos2 + slack >= targ2)
    
    gene_df2 <- glazer_genes %>%
      filter(chr == chr_target) %>%
      filter(pos1 - slack <= targ2, pos2 + slack >= targ2)
    
    gene <- unique(gene_df1$gene_id, gene_df2$gene_id)
    gene_df <- data.frame(chr = chr_target, window = pos_target, gene)
    
  } else{
    
    chr_target <- sites_df$chr[index]
    pos_target <- sites_df$pos[index]
    
    gene_df <- glazer_genes %>%
      filter(chr == chr_target) %>%
      filter(pos1 - slack <= pos_target, pos2 + slack >= pos_target)
    
    
    gene <- unique(gene_df$gene_id)
    gene_df <- data.frame(chr = chr_target, pos = pos_target, gene)
    
  }
  
}

tmp_wind <- sites_df %>% filter(chr == 1) %>% select(window) %>% unlist


lapply()
find_gene_by_chr_pos(1, tmp_wind[50], glazer_genes)
