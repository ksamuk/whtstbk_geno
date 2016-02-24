install.packages("PopGenome")
library("PopGenome")

wht_vcf <- readData(path = "data/vcf/wht.vcf", format = "VCF")

wht_vcf <- readVCF("data/vcf/whtstbk_nova_scotia_filtered.vcf.gz", frompos = 1, topos = 8000000, tid = "chrV", numcols = 10000)

wht_vcf <- neutrality.stats(wht_vcf)
get.neutrality(wht_vcf)[[1]]
