#### Subsets, filters and reformats a VCF 
#### end point is a 'geno_df' 
#### example run: Rscript 01_process_vcf_to_geno_df.R data/vcf/whtstbk_master.vcf.bgz multi 2

################################################################################
# Check and parse command arguments (if any)
################################################################################

if (length(commandArgs(TRUE)) < 0){
  
  stop("Invalid arugments! Run like Rscript 01_process_vcf_to_geno_df.R data/vcf/vcf.file.vcf multi|single [cores]")
  
}

args <- commandArgs(TRUE)
args <- as.list(args)

vcf_file <- args[[1]]
vcf_index <- paste0(vcf_file, ".tbi")

paral <- args[[2]]

################################################################################
# Libraries
################################################################################

# source("http://bioconductor.org/biocLite.R")
# biocLite("BiocUpgrade")
# source("http://bioconductor.org/workflows.R")
# workflowInstall("variants")

library(VariantAnnotation)
library(dplyr)
library(tidyr)
library(ggplot2)
library(cowplot)
library(ggthemes)
library(magrittr)

list.files("functions", full.names = TRUE) %>% sapply(.,source, verbose = FALSE, echo = FALSE) %>% invisible

################################################################################
# Main processing
################################################################################

#index the vcf
indexTabix(file = vcf_file, format = "vcf")

# list of contigs to process
chr_list <- c(paste0("chr", as.roman(1:21)), "chrUn")

# apply geno_df function to vcf file
# parses args for multi/single core action

if (paral == "single"){
  lapply(chr_list, process_contig_to_geno_df, vcf_file = vcf_file, vcf_index = vcf_index, overwrite = FALSE)
} else{
  mclapply(chr_list, process_contig_to_geno_df, 
           vcf_file = vcf_file, vcf_index = vcf_index, overwrite = FALSE, 
           mc.cores = args[[3]], mc.silent = FALSE, mc.preschedule = FALSE)
}




