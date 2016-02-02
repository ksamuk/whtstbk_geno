read_structure_pop <- function(str_file_name, meta_df){
  
  # read in structure file and set column names
  pop_name <- gsub("[^A-Z]|meanQ","", str_file_name)
  df <- read.table(str_file_name)
  max_k <- as.numeric(gsub("\\.[0-9]\\.","", str_file_name) %>% gsub("[^0-9]", "", .))
  k_cols <- paste0("k",1:max_k)
  names(df) <- k_cols
  
  # join in str_df (scoped)
  meta_df <- meta_df %>% filter(pop == pop_name)
  df <- data.frame(meta_df, df)
  
  # clean ids
  df$id <- df$id %>% 
    gsub("whtstbk_gbs_2012_brds_", "2012_", .) %>%
    gsub("\\.", "_", .) %>%
    str_trim
  
  #add 
  
  df %>% 
    gather_(key_col = "k", value_col = "q.value", 
            gather_cols = k_cols) %>%
    arrange(id, k)

}