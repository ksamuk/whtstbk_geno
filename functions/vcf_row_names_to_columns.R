# splits a vcf row names into four separate columns

vcf_row_names_to_columns <- function(x){
	
	x_split <- strsplit(x, split = ":|_|/") 
	
}