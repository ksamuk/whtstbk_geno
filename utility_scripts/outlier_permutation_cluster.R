# outlier window permutation test

library("parallel")

is.outlier <- function(x, cutoff = 0.95){
  
  x95 <- quantile(x, na.rm = TRUE, probs = cutoff)[1]
  return(x >= x95)
}


#fx_df <- read.table("data/stats/snp_stats_master.txt", header = TRUE, stringsAsFactors = FALSE)

fx_df <- read.table("snp_stats_master.txt", header = TRUE, stringsAsFactors = FALSE)

outlier_wht_cmn <- is.outlier(fx_df$fst_wht_cmn, cutoff = 0.99)
outlier_wht_cmn2012 <- is.outlier(fx_df$fst_wht_cmn2012, cutoff = 0.99)
outlier_wht_cbr <- is.outlier(fx_df$fst_wht_cbr, cutoff = 0.99)
outlier_cbr_cmn <- is.outlier(fx_df$fst_cbr_cmn, cutoff = 0.99)

#fx_df$outlier_fst_joint <- (outlier_wht_cmn | outlier_wht_cbr) & !(outlier_cbr_cmn)
#fx_df$outlier_xtx_joint <- with(fx_df, (outlier_xtx_wht_cbr | outlier_xtx_wht_cmn) & !(outlier_xtx_cmn_cbr))

# extract each outlier stat into its own outlier df (in a list)

# they say the classis never go out of style but...they do. they do.

outlier_names <- c("fst_wht_cmn", "fst_wht_cbr", "fst_cbr_cmn", "xtx_wht_cbr", "xtx_wht_cmn", "xtx_cmn_cbr", "fst_wht_cmn2012", "xtx_wht_cmn2012")

outlier_df_list <- list()

outlier_df_list[[1]]  <- fx_df %>%
  select(chr, pos)%>%
  mutate(outlier = outlier_wht_cmn) %>%
  filter(!is.na(outlier)) %>%
  arrange(chr, pos)

outlier_df_list[[2]]  <- fx_df %>%
  select(chr, pos)%>%
  mutate(outlier = outlier_wht_cbr) %>%
  filter(!is.na(outlier)) %>%
  arrange(chr, pos)

outlier_df_list[[3]]  <- fx_df %>%
  select(chr, pos)%>%
  mutate(outlier = outlier_cbr_cmn) %>%
  filter(!is.na(outlier)) %>%
  arrange(chr, pos)

outlier_df_list[[4]]  <- fx_df %>%
  select(chr, pos, outlier_xtx_wht_cbr)%>%
  rename(outlier = outlier_xtx_wht_cbr) %>%
  filter(!is.na(outlier)) %>%
  arrange(chr, pos)

outlier_df_list[[5]]  <- fx_df %>%
  select(chr, pos, outlier_xtx_wht_cmn)%>%
  rename(outlier = outlier_xtx_wht_cmn) %>%
  filter(!is.na(outlier)) %>%
  arrange(chr, pos)

outlier_df_list[[6]]  <- fx_df %>%
  select(chr, pos, outlier_xtx_cmn_cbr)%>%
  rename(outlier = outlier_xtx_cmn_cbr) %>%
  filter(!is.na(outlier)) %>%
  arrange(chr, pos)

outlier_df_list[[7]]  <- fx_df %>%
  select(chr, pos)%>%
  mutate(outlier = outlier_wht_cmn2012) %>%
  filter(!is.na(outlier)) %>%
  arrange(chr, pos)

outlier_df_list[[8]]  <- fx_df %>%
  select(chr, pos, outlier_xtx_wht_cmn.2012)%>%
  rename(outlier = outlier_xtx_wht_cmn.2012) %>%
  filter(!is.na(outlier)) %>%
  arrange(chr, pos)

rm(fx_df)

# assign windows

permute_outlier_windows <- function(index, outlier_df_list, outlier_names = outlier_names, reps = 10000, cutoff = 0.95, fast = TRUE){
  
  print(paste0("Processing ", outlier_names[index], "..."))
  
  outlier_df <- outlier_df_list[[index]]
  
  # create windows  
  outlier_df$w_pos2 <- (((outlier_df$pos / 75000) %>% floor) + 1)*75000
  outlier_df$window_start <- outlier_df$w_pos2 - 74999
  
  outlier_df <- outlier_df %>%
    select(chr, window_start, pos, outlier) %>%
    mutate(chr_window = paste0(chr, "_", window_start))
  
  # toss windows with fewer than 3 snps
  
  windows_to_keep <- outlier_df %>%
    group_by(chr_window) %>%
    tally %>%
    filter(n >= 3) %>%
    select(chr_window) %>% unlist
  
  outlier_df <- outlier_df %>%
    filter(chr_window %in% windows_to_keep)

  # perform permutation
  
  permute_outlier <- . %>% 
    mutate(outlier = sample(outlier, replace = FALSE)) %>%
    group_by(chr_window) %>%
    summarise(outlier_count = sum(outlier))
  
  perm_df <- replicate(reps, permute_outlier(outlier_df), simplify = FALSE)
  perm_df <- rbind_all(perm_df)
  
  if (fast == TRUE){
    
    cutoff_df <- perm_df %>% 
      arrange(chr_window) %>%
      group_by(chr_window) %>%
      summarise(cutoff = quantile(rnorm(outlier_count), probs =  cutoff) %>% as.numeric())
    
    emp_df <- outlier_df %>%
      group_by(chr_window) %>%
      summarise(outlier_count = sum(outlier)) %>%
      left_join(cutoff_df)
    
    emp_df$outlier_window <- emp_df$outlier_count >= emp_df$cutoff
    
    names(emp_df)[4] <- outlier_names[index]
    
    emp_df[,c(1,4)]
    
  } else{
    
    emp_df <- outlier_df %>%
      group_by(chr_window) %>%
      summarise(outlier_count = sum(outlier))
    
    calc_permutation_pval <- function(window, emp_df, perm_df){
      
      perm_vec <- perm_df %>%
        filter(chr_window ==  window) %>%
        select(outlier_count) %>%
        unlist %>% as.numeric
      
      emp_count <- emp_df %>%
        filter(chr_window ==  window) %>%
        select(outlier_count) %>%
        unlist %>% as.numeric
      
      (sum((perm_vec >= emp_count))+ 1)/(length(perm_vec) +1)
      
    }
    
    pval <- lapply(emp_df$chr_window, calc_permutation_pval, emp_df, perm_df) %>% unlist
    
    emp_df$pval <- pval
    names(emp_df)[3] <- paste0(outlier_names[index], "_pval")
    
    names(emp_df)[2] <- paste0(outlier_names[index], "_count")
    
    emp_df
    
  }
  
  
}

indexes <- c(1:8)

# outlier_df <- mclapply(indexes, permute_outlier_windows, reps = 10000, cutoff = 0.99, outlier_names = outlier_names, outlier_df_list = outlier_df_list, fast = FALSE,
#                        mc.preschedule = FALSE, mc.set.seed = TRUE,
#                        mc.silent = FALSE, mc.cores = 6,
#                        mc.cleanup = TRUE, mc.allow.recursive = TRUE)

#outlier_df <- lapply(indexes, permute_outlier_windows, reps = 10, cutoff = 0.99, outlier_names = outlier_names, outlier_df_list = outlier_df_list, fast = FALSE)

outlier_df <- mclapply(indexes, permute_outlier_windows, reps = 10000, cutoff = 0.99, outlier_names = outlier_names, outlier_df_list = outlier_df_list, fast = FALSE, mc.preschedule = FALSE, mc.set.seed = TRUE, mc.silent = FALSE, mc.cores = 8, mc.cleanup = TRUE, mc.allow.recursive = TRUE)

outlier_out_df <- outlier_df[[1]] %>% left_join(outlier_df[[2]], by = c("chr_window")) %>%
  left_join(outlier_df[[3]], by = c("chr_window")) %>%
  left_join(outlier_df[[4]], by = c("chr_window")) %>%
  left_join(outlier_df[[5]], by = c("chr_window")) %>%
  left_join(outlier_df[[6]], by = c("chr_window")) %>%
  left_join(outlier_df[[7]], by = c("chr_window")) %>%
  left_join(outlier_df[[8]], by = c("chr_window")) %>%
  arrange(chr_window)

chr <- strsplit(outlier_out_df$chr_window, split = "_") %>% lapply(function(x)x[[1]]) %>% unlist %>% as.numeric
pos <- strsplit(outlier_out_df$chr_window, split = "_") %>% lapply(function(x)x[[2]]) %>% unlist %>% as.numeric

outlier_out_df <- data.frame(chr, pos, outlier_out_df[,-c(1)]) %>%
  arrange(chr, pos)

write.table(outlier_out_df, "outlier_permutation_df.txt", quote = FALSE, row.names = FALSE)