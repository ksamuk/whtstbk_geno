# Libraries ---------------------------------------------------------------

source("http://bioconductor.org/biocLite.R")
biocLite("BiocUpgrade")
source("http://bioconductor.org/workflows.R")
workflowInstall("variants")

library(VariantAnnotation)
library(dplyr)

# read in a vcf file
vcf_dat <- readVcf("data/test.vcf", genome = "glazer")
vcf_index <- TabixFile(file = "data/test.vcf.idx")


# compression and indexing -------------------------------------------------

#bgzip a vcf
bgzip("data/test.vcf")

#index a bgziped vcf
indexTabix(file = "data/test.vcf.bgz", format = "vcf")
tab <- TabixFile("data/test.vcf.bgz", "data/test.vcf.bgz.tbi")

# extract a subset of a vcf -----------------------------------------------
vcf <- readVcf(tab, "glazer")

# specific chromosome/contig

chrI <- vcf[seqnames(rowRanges(vcf)) == "chrI"]
geno(chrI)
vcfFixed(chrI)

# subset of individuals

# the header (duh)
header(vcf_dat)

# the list of samples ooo shiny
samples(header(vcf_dat))

# sort of like a tbl_df header
info(vcf_dat)

# genotype extraction?
geno(vcf_dat)$GT

gt_dat <- readGT("data/test.vcf", nucleotides = TRUE)