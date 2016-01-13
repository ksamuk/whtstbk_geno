fst_per_locus <- function(fstat_df, pop1, pop2){
  
  # filter fstat data for target populations
  df_sub <- fstat_df %>%
    filter(pop == pop1 | pop == pop2)%>%
    droplevels
  
  # convert fstat_df column names to chr/pos columns
  
  perloc_df <- fstat_label_to_columns(df_sub[,-1])
  perloc_df$fst <- wc(df_sub)$per.loc$FST
  perloc_df <- data.frame(comparison = paste0(pop1, "_", pop2), perloc_df)
  perloc_df
  
}
