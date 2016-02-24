read_structure <- function(str_file_name, ids, meta_df){
  
  # read in structure file and set column names
  str_df <- read.table(str_file_name)
  max_k <- as.numeric(gsub("[^0-9]*|2014","", str_file_name)) %>% substr(1,1)
  k_cols <- paste0("k", 1:max_k)
  names(str_df) <- k_cols
  str_df <- data.frame(id = unique(ids[,1]), str_df)
  
  # join in str_df (scoped)
  str_df <- left_join(str_df, meta_df)
  
  # clean ids
  str_df$id <- str_df$id %>% 
    gsub("whtstbk_gbs_2012_brds_", "2012_", .) %>%
    gsub("\\.", "_", .) %>%
    str_trim
  
  #add 
  
  str_df %>% 
    gather_(key_col = "k", value_col = "q.value", 
            gather_cols = k_cols)
}