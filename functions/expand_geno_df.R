# expands / injects extra information into a geno_df

expand_geno_df <- function(geno_df){
  
  # create separate genotype fields
  geno_df$geno_1 <- strsplit(geno_df$GT, split = "/") %>% lapply(function(x)x[1]) %>% as.numeric
  geno_df$geno_2 <- strsplit(geno_df$GT, split = "/") %>% lapply(function(x)x[2]) %>% as.numeric
  
  # coerce ref/alt to character
  geno_df$alt <- geno_df$alt %>% as.character
  geno_df$ref <- geno_df$ref %>% as.character
  
  # map to the actual allelic identities
  
  # split the alt column
  geno_df$alt_split <- geno_df$alt %>% unlist %>% strsplit(split = ",")
  
  print("converting genotypes to nucleotides...")
  # convert genotypes to nucleotides
  geno_df <- geno_df %>%
    rowwise %>%
    mutate(geno_1_nuc = base::ifelse(geno_1 == 0, ref,  alt_split[geno_1])) %>% 
    mutate(geno_2_nuc = base::ifelse(geno_2 == 0, ref,  alt_split[geno_2]))
  
  # split pop from ID
  geno_df$pop <- base::gsub("[^A-Z]", "", geno_df$id)
  geno_df$ind <- base::gsub("[^[0-9]*|2012|2014", "", geno_df$id)
 
  # year
  geno_df$year <- ifelse(grepl("whtstbk", as.character(geno_df$id)), 2012, 2014)
                         
  geno_df <- with(geno_df, data.frame(year, pop, ind,chrom, pos, ref, alt, GT, geno_1, geno_2, geno_1_nuc, geno_2_nuc))
  
  geno_df
  
} 