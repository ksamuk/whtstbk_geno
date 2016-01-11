#### Subsets, filters and reformats a VCF 
#### end point is a 'geno_df' 

process_contig_to_geno_df <- function(vcf_file = vcf_file, overwrite = FALSE, chr){
  
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
  
  vcf_sub <- read_vcf_subset(vcf_file, vcf_index, chr = chr)
  
  # filter vcf (site metrics) ------------------------------------------------------------
  
  # read in info fields
  vcf_info <- vcf_info_df(vcf_sub)
  
  # plot qc summaries
  # produces ggplot warnings, safe to ignore them
  qc_figure <- vcf_info %>%
    plot_qc_figures(frac = 1)
  
  ggsave(paste0("reports/", slug, "/", "qc_plot_", chr, "_pre.png"), plot = qc_figure, width = 11, height = 8.5)
  
  # apply a QD and MQ filter
  # also plot results
  vcf_info <- vcf_info %>% 
    filter(QD >= 20, MQ >= 40)
  
  qc_figure <- vcf_info %>%
    plot_qc_figures(frac = 1)
  
  ggsave(paste0("reports/", slug, "/", "qc_plot_", chr, "_post.png"), plot = qc_figure, width = 11, height = 8.5)
  
  # apply site filter to genotypes -----------------------------------------------
  
  # extract genotype fields and format as dataframes
  vcf_geno <- suppressWarnings(vcf_geno_df(vcf_sub))
  
  # apply site filter and:
  # - remove NAs (and '.'s)
  # - depth >= 8
  # - genotype quality >= 30
  
  vcf_geno <- vcf_geno %>%
    filter(pos %in% vcf_geno$pos) %>%
    filter(GT != ".") %>%
    filter(DP >= 8) %>%
    filter(GQ >= 30)
  
  # split GT field into two columns
  
  vcf_geno <- expand_geno_df(vcf_geno)
  
  # write to file ---------------------------------------------------------
  
  write.table(vcf_geno, file = gzfile(paste0("data/geno_df/", chr, ".gz")))
  }
}