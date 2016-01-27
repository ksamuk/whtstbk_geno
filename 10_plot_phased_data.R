#### Subsets, filters and reformats a VCF 
#### end point is a 'geno_df' 
#### example run: Rscript 01_process_vcf_to_geno_df.R data/vcf/whtstbk_master.vcf.bgz multi 2

################################################################################
# Libraries
################################################################################

list.files("functions", full.names = TRUE) %>% sapply(.,source, verbose = FALSE, echo = FALSE) %>% invisible

file_list <- list.files("data/phased_vcf", pattern=".vcf$", full.names = TRUE)

vcf <- file_list[1]

tmp <- phased_vcf_to_geno_df(file_list[1])
