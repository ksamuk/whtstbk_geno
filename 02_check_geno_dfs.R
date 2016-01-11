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
library(readr)

list.files("functions", full.names = TRUE) %>% sapply(.,source, verbose = FALSE, echo = FALSE) %>% invisible

# read in a geno_df

geno_dfs <- list.files("data/geno_df", full.names = TRUE)

gdf <- read.table(file = geno_dfs[10], stringsAsFactors = FALSE, header = TRUE)
