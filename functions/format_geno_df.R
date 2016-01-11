format_geno_df <- function(geno_matrix, name = "value"){
	
	# covert rownames to columns
	geno_df <- geno_matrix %>% 
		rownames %>% 
	  vcf_row_names_to_columns %>%
		#lapply(vcf_row_names_to_columns) %>% 
		do.call("rbind", .) %>%
	  data.frame
	
	names(geno_df) <- c("chrom", "pos", "ref", "alt")
	
	# gather geno fields into columns
	geno_df <- gather(data.frame(geno_df, geno_matrix), id, geno, -chrom, -pos, -ref, -alt)
	names(geno_df)[6] <- name
	geno_df
}
