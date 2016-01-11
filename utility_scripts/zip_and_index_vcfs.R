#### (compress) and index a folder of vcf files 

################################################################################
# Libraries
################################################################################

# source("http://bioconductor.org/biocLite.R")
# biocLite("BiocUpgrade")
# source("http://bioconductor.org/workflows.R")
# workflowInstall("variants")

library(VariantAnnotation)
library(dplyr)
library(magrittr)

list.files("functions", full.names = TRUE) %>% sapply(.,source, verbose = FALSE, echo = FALSE) %>% invisible

################################################################################
# Main processing
################################################################################

# the list of vcfs to compress / index
file_list <- list.files("data/vcf/by_chromo", pattern = ".gz$",full.names = TRUE)

#index the vcf

for (i in file_list){
  #zipped_vcf <- bgzip(i, overwrite = TRUE)
  indexTabix(file = i, format = "vcf")
}
