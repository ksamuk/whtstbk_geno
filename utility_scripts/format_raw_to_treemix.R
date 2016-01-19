# format a .raw plink file to treemix format
# assumes plink is malformed from vcf conversion
# e.g. the FID and IID columns are identical

#library("data.table")

library("readr")
library("stringr")
library("dplyr")
library("tidyr")

list.files("functions", full.names = TRUE) %>% sapply(.,source, verbose = FALSE, echo = FALSE) %>% invisible

plink_file <- list.files("data/treemix", pattern = ".raw", full.names = TRUE)
raw_file <- read.table(plink_file, header = TRUE, stringsAsFactors = FALSE)
meta_df <- read.csv("metadata/mega_meta.csv")

# fix population names
raw_file$FID <- gsub("[^A-Z]*", "", raw_file$FID)

# filter out low sample size pops
pop_counts <- meta_df %>% group_by(pop) %>% tally %>% filter(n > 4)
meta_df <- meta_df %>% 
  filter(pop %in% pop_counts$pop)

raw_file <- raw_file %>%
  filter(FID %in% pop_counts$pop)

# collapse CP pops
raw_file$FID[raw_file$FID == "CPN"|raw_file$FID == "CPSE"] <- "CP"

# add in species labels
pop_suffixes <- left_join(data.frame(id = raw_file$IID), meta_df)
raw_file$FID <- paste0(raw_file$FID, "_", meta_df$species)

# filter out 2012 data
raw_file <- raw_file %>% 
  filter(!grepl("2012", IID))

raw_file <- raw_file %>%
  select(-PAT, -MAT, -SEX, -PHENOTYPE, -IID)

# remove pops/individuals

# remove 2012 individuals
#  raw_file <- raw_file %>%
#  filter(!grepl("2012", IID)) 

# tally allele counts in populations

# convert to long
raw_long <- gather(raw_file, key = locus, value = allele, -FID)

# recover memory
rm(raw_file)

pop_sums <- raw_long %>%
  group_by(FID, locus, allele) %>%
  tally

# recover memory
rm(raw_long)

# split chr labels into chr pos + convert to numeric
chr_pos <- pedraw_label_to_columns(as.character(pop_sums$locus))
chr_pos$chr <- as.character(chr_pos$chr)
chr_pos$chr[chr_pos$chr == "chrUn"] <- "chrXXII"

chr_pos$chr <- chr_pos$chr %>% gsub("chr", "", .) %>% as.roman %>% as.numeric

# re gather data frame
pop_sums <- data.frame(pop = pop_sums$FID, chr_pos, allele = pop_sums$allele, n = pop_sums$n)

# recover memory
rm(chr_pos)

# arrange by chr/pos
pop_sums <- pop_sums %>%
  arrange(chr, pos)

# format:
# Han  Sardinian French Karitiana Yoruba 
# 2,66 2,54      4,52   0,24      21,21

# generate the allele count strings
pop_count <- pop_sums %>%
  group_by(pop, chr, pos) %>%
  summarise(count_out = paste0(sum(n[allele == 0]*2, n[allele == 1], na.rm = TRUE), "," , sum(n[allele == 2]*2, n[allele == 1], na.rm = TRUE)))

# coerce pos back to numeric
pop_count$pos <- as.numeric(as.character(pop_count$pos))

# spread the data and remove chr pos
pop_count <- spread(pop_count, pop, count_out) %>%
  arrange(chr, pos) %>%
  select(-chr, -pos)

# out file names
out_file <- gsub(".raw", ".treemix.txt", plink_file)

# write to file
write.table(pop_count, file = out_file, row.names = FALSE, quote = FALSE)

