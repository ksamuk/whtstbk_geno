#### Subsets, filters and reformats a VCF 
#### end point is a 'geno_df' 
#### example run: Rscript 01_process_vcf_to_geno_df.R data/vcf/whtstbk_master.vcf.bgz multi 2

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
library(stringr)

list.files("functions", full.names = TRUE) %>% sapply(.,source, verbose = FALSE, echo = FALSE) %>% invisible

################################################################################
# Main processing
################################################################################

#list of files
file_list <- list.files("data/vcf/by_chromo", pattern=".gz$", full.names = TRUE)

# apply geno_df function to vcf file
# parses args for multi/single core action

lapply(file_list, process_contig_to_geno_df_by_chromo, overwrite = FALSE)
