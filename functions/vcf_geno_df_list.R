# return a formatted list of genotype fields
# ('long' dataframe)

vcf_geno_df <- function(vcf){
	
	# get geno field matrices
	vcf_geno <- geno(vcf) %>% as.list
	
	# format the GT matrix
	# includes step of splitting names
	geno_names_gt <- format_geno_df(vcf_geno[[1]], name = "GT")
	
	# format remaining matrices (no names)
	geno_names <- names(vcf_geno)[c(2,3)]
	vcf_geno <- vcf_geno[c(2,3)]
	
	# format as dataframes 
	df_list <- Map(format_geno_df_no_rows, vcf_geno, names(vcf_geno))
	
	geno_df <- data.frame(geno_names_gt, df_list[[1]][2], df_list[[2]][2])
	geno_df$pos <- geno_df$pos %>% as.character %>% as.numeric
	geno_df
	
}