#### Subsets, filters and reformats a VCF 
#### end point is a 'geno_df' 

process_contig_to_geno_df_by_chromo <- function(vcf_file, overwrite = FALSE){
  
  #set fake chr for compatibility
  chr_position <- gregexpr("chr[IXVUn]+", vcf_file) %>% unlist
  chr_length <- attributes(gregexpr("chr[IXVUn]+", vcf_file)[[1]])$match.length
  chr <- base::substr(vcf_file, chr_position, chr_position+(chr_length-1))
  
  # initialize index file and output file name
  vcf_index <- paste0(vcf_file, ".tbi")
  out_file_name <- paste0("data/geno_df/", chr, ".gz")
  
  # unqie slug for filenames
  slug <- vcf_file %>% 
    strsplit(split = "/") %>% 
    lapply(function(x)x[3]) %>% 
    unlist %>% 
    gsub(".vcf.bgz", "", .)
  
  dir.create(paste0("reports/", slug))
  
  if(overwrite == FALSE & file.exists(out_file_name)){
    
    print(paste0(out_file_name, " exists, skipping..."))
    
  } else{
  
  print(paste0("processing ", chr, "..."))
  
  # extract a subset of a vcf -----------------------------------------------
  
  param <- ScanVcfParam(geno = c("GT", "DP", "GQ"), info = c("DP", "MQ", "QD", "AN", "InbreedingCoeff"))
  print("reading vcf...")
  tab <- TabixFile(vcf_file, vcf_index)  
  vcf_sub <- readVcf(tab, genome = "", param)
  
  # filter vcf (site metrics) ------------------------------------------------------------
  print("harvesting info fields...")
  # read in info fields
  vcf_info <- vcf_info_df(vcf_sub)
  
  print("performing PRE qc...")
  # plot qc summaries
  # produces ggplot warnings, safe to ignore them
  qc_figure <- vcf_info %>%
    plot_qc_figures(frac = 1)
  
  ggsave(paste0("reports/", slug, "/", "qc_plot_", chr, "_pre.png"), plot = qc_figure, width = 11, height = 8.5)
  
  print("harvesting ALT field...")
  #harvest the ALT field (ensures polyallelic SNPs are propagated)
  vcf_info$alt <- mcols(vcf_sub)$ALT %>% 
    lapply(str_c, collapse = ",")
  
  # apply a QD and MQ filter
  # also plot results
  vcf_info <- vcf_info %>% 
    mutate(PI = AN/max(AN)) %>% 
    filter(QD >= 20, MQ >= 40, InbreedingCoeff > -1, PI >= 0.50)
  
  print("performing POST qc...")
  
  qc_figure <- vcf_info %>%
    plot_qc_figures(frac = 1)
  
  ggsave(paste0("reports/", slug, "/", "qc_plot_", chr, "_post.png"), plot = qc_figure, width = 11, height = 8.5)
  
  # apply site filter to genotypes -----------------------------------------------
  
  print("extracting genotype fields...")
  # extract genotype fields and format as dataframes
  vcf_geno <- suppressWarnings(vcf_geno_df(vcf_sub))
  rm(vcf_sub)
  
  # apply site filter and
  vcf_geno <- vcf_geno %>%
    filter(pos %in% vcf_info$pos)
  
  # add in long form alt field
  vcf_geno$alt <- vcf_info$alt
  
  # genotype filters:
  # - remove NAs (and '.'s)
  # - depth >= 8
  # - genotype quality >= 30
  print("filtering genotypes...")
  vcf_geno <- vcf_geno %>% 
    filter(GT != ".") %>%
    filter(DP >= 8) %>%
    filter(GQ >= 30) %>%
    select(-DP, -GQ)
  
  # split GT field into two columns
  print("expanding geno_df...")
  vcf_geno <- expand_geno_df(vcf_geno)
  vcf_geno$pop <- as.factor(vcf_geno$pop)
  vcf_geno$year <- as.factor(vcf_geno$year)
  
  # write to file ---------------------------------------------------------
  
  print("writing to file...")
  write.table(vcf_geno, file = gzfile(paste0("data/geno_df/", chr, ".gz")))
  rm(vcf_geno)
  rm(vcf_info)
  }
}