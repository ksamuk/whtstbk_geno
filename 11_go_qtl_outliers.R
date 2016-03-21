# GO / QTL Anaylsis

################################################################################
# initials
################################################################################

library("dplyr")
<<<<<<< HEAD
library("ggplot2")
library("ggrepel")
library("biomaRt")

################################################################################
# raw data
################################################################################

library("parallel")

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

glazer_genes <- glazer_genes %>%
  group_by(chr, gene_id) %>%
  summarize(pos1 = min(pos1), pos2 = max(pos2)) %>%
  ungroup

################################################################################
# joining data
################################################################################



find_gene_by_chr_pos <- function(index, glazer_genes, slack = 0, window_size = 75000){
  
  if(window_size > 0){
    
    chr_target <- sites_df$chr[index]
    pos_target <- sites_df$window[index]
    
    targ1 <- pos_target
    targ2 <- pos_target + (window_size - 2)
    
    sites <- targ1 : targ2
    
    gene_df1 <- glazer_genes %>%
      filter(chr == chr_target) %>%
      filter(pos1 > targ1) %>%
      filter(pos2 < targ2) %>%
      group_by(gene_id) %>%
      summarise(overlap = (any(sites %in% pos1:pos2))) %>%
      filter(overlap == TRUE)
      #filter(pos1 - slack <= targ1, pos2 + slack >= targ2)
    
    if(nrow(gene_df1) > 0){
      data.frame(chr = chr_target, pos1 = targ1, pos2 = targ2, gene = unique(gene_df1$gene_id))
    }else{
      data.frame(chr = chr_target, pos1 = targ1, pos2 = targ2, gene = NA)
    }
    
    
  } else{
    
    chr_target <- sites_df$chr[index]
    pos_target <- sites_df$pos[index]
    
    gene_df <- glazer_genes %>%
      filter(chr == chr_target) %>%
      filter(pos1 - slack <= pos_target, pos2 + slack >= pos_target)
    
    gene <- unique(gene_df$gene_id)
    if(length(gene) == 0){
      gene <- NA
    }
    gene_df <- data.frame(chr = chr_target, pos = pos_target, gene)
    gene_df
    
  }
  
}

match_df <- list()

for (i in 1:length(unique(wind_df$chr))){
  
  chr_target <- unique(wind_df$chr)[i]
  
  cat(paste0("chr ", chr_target, "..."))
  
  windows <- wind_df %>%
    filter(chr == chr_target) %>%
    select(window) %>%
    unlist %>% as.numeric
  
  match_df[[i]] <- lapply(windows, find_gene_by_chr_pos, chr_target = chr_target, glazer_genes = glazer_genes)
}

match_df <- bind_rows(unlist(match_df, recursive = FALSE))


tmp <- setNames(glazer_genes, c("chr", "gene", "gene_pos1", "gene_pos2"))
match_df <- left_join(match_df, tmp, by = c("chr", "gene"))
tmp <- setNames(wind_df, c("chr", "pos1", "fst_outlier_adj", "xtx_outlier_adj"))
match_df <- left_join(match_df, tmp, by = c("chr", "pos1"))

################################################################################
# finding annotation data using biomaRt
################################################################################

# define biomart object for gacu
ensembl_mart <- useMart("ensembl")
gacu_mart <- useMart("ensembl",dataset="gaculeatus_gene_ensembl")

# full lists of filters/attributes
attributes <-  listAttributes(gacu_mart, what=c("name","description","page"))
filters <-  listFilters(gacu_mart, what=c("name","description"))

#View(attributes)
#View(filters)

# first, get gene IDs (glazer provides two transcript ids instead of gene + transcript ids)

gene_df <- getBM(attributes = c("ensembl_gene_id", 
                                 "ensembl_transcript_id"), mart = gacu_mart, filters = "ensembl_transcript_id", values = match_df$gene)
names(gene_df) <- c("gene", "transcript")
names(match_df)[4] <- "transcript"
match_df <- left_join(match_df, gene_df)

# get ensembl annotation info
annot_df <- getBM(attributes = c("ensembl_gene_id", 
                                 "ensembl_transcript_id",
                                 "external_gene_name",
                                 "go_id",
                                 "name_1006",
                                 "definition_1006",
                                 "go_linkage_type",
                                 "namespace_1003",
                                 "goslim_goa_accession",
                                 "goslim_goa_description"), 
                  mart = gacu_mart, filters = "ensembl_gene_id", values = match_df$gene)

write.table(annot_df, "data/genes/gacu_annotations.txt", quote = FALSE, row.names = FALSE)
write.table(match_df, "data/genes/outlier_genes.txt", quote = FALSE, row.names = FALSE)
  
################################################################################
# plots
################################################################################

match_df %>% 
  filter(chr == 1) %>%
  filter(fst_outlier_adj == TRUE) %>%
  ggplot(aes(x = pos1, xmin = pos1, xmax = pos2, ymin = 1, ymax = 2, fill = fst_outlier_adj))+
  geom_rect()+
  geom_rect(aes(x = gene_pos1, xmin = gene_pos1, xmax = gene_pos2, ymin = 2, ymax = 3))+
  #geom_text(aes(x = gene_pos1, label = gene, y = 3, check_overlap = TRUE), size = 2)+
  geom_label_repel(aes(x = gene_pos1, label = gene, y = 3), size = 2)+
  scale_fill_manual(values = c("grey", "red"))
