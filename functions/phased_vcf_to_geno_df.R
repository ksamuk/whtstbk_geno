# return a formatted list of genotype fields
# ('long' dataframe)

phased_vcf_to_geno_df <- function(vcf){
  
  bgzip(vcf, overwrite = TRUE)
  vcf <- paste0(vcf, ".bgz")
  tab <- TabixFile(vcf, paste0(vcf, ".tbi"))  
  vcf <- readVcf(tab, genome = "", param)
	
	# get geno field matrices
	vcf_geno <- geno(vcf) %>% as.list
	
	# format the GT matrix
	# includes step of splitting names
	geno_df <- format_geno_df(vcf_geno[[1]], name = "GT")
	
	geno_df$pos <- geno_df$pos %>% as.character %>% as.numeric
	geno_df
	
}