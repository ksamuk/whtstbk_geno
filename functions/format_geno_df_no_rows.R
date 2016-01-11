format_geno_df_no_rows <- function(geno_matrix, name = "value"){
	
	# gather geno fields into columns
  geno_matrix <- data.frame(geno_matrix)
	geno_df <- gather(geno_matrix, id, geno)
	names(geno_df)[2] <- name
	
	geno_df
}
