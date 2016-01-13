read_structure <- function(str_file_name, meta_df){
  
  # read in structure file and set column names
  df <- read.table(str_file_name)
  max_k <- as.numeric(gsub("[^0-9]*","", str_file_name))
  k_cols <- paste0("k",1:max_k)
  names(df) <- k_cols
  
  # join in str_df (scoped)
  df <- data.frame(meta_df, df)
  
  # clean ids
  df$id <- df$id %>% 
    gsub("whtstbk_gbs_2012_brds_", "2012_", .) %>%
    gsub("\\.", "_", .) %>%
    str_trim
  
  #add 
  
  df %>% 
    gather_(key_col = "k", value_col = "q.value", 
            gather_cols = k_cols)
}