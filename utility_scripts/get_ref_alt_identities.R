# harvest allelic identities from a VCF file

#install.packages("vcfR")
library("vcfR")
library("dplyr")

vcf <- read.vcfR("data/vcf/whtstbk_bial_nomaf_nosex.vcf.gz")

# the fixed fields (chr pos ref alt...)

fix_df <- vcf@fix
rm(vcf)

fix_df <- fix_df %>%
  data.frame %>%
  select(CHROM, POS, REF, ALT) 

names(fix_df) <- c("chr", "pos", "ref", "alt")

write.table(fix_df, "metadata/ref_alt_nova_scotia_master.txt", quote = FALSE, row.names = FALSE)
