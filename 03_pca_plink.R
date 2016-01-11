# plot PCA of genotype data

library(adegenet)
library(pegas)

plink_file <- list.files("data/other_formats", pattern = ".raw$", full.names = TRUE)

wht_gen <- read.PLINK(plink_file, parallel = FALSE)

