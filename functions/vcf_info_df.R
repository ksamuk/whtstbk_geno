# return a formatted list of site info fields
# ('long' dataframe)

vcf_info_df <- function(vcf){
	
	# split row names in to columns
	vcf_info <- info(vcf) %>% 
		rownames %>% 
		vcf_row_names_to_columns %>% 
		do.call("rbind", .) %>%
	  data.frame
	
	names(vcf_info) <- c("chrom", "pos", "ref", "alt")
	
	# return joined df
	vcf_info <- data.frame(vcf_info, info(vcf))
	vcf_info$pos <- as.numeric(as.character(vcf_info$pos))
	vcf_info
}


