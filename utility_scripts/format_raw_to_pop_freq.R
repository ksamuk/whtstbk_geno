# format a .raw plink file to treemix format
# assumes plink is malformed from vcf conversion
# e.g. the FID and IID columns are identical

#library("data.table")

library("readr")
library("stringr")
library("dplyr")
library("tidyr")

list.files("functions", full.names = TRUE) %>% sapply(.,source, verbose = FALSE, echo = FALSE) %>% invisible

plink_file <- list.files("data/other_formats", pattern = "full", full.names = TRUE)
#plink_file <- list.files("data/treemix", pattern = ".raw", full.names = TRUE)
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
# raw_file <- raw_file %>% 
#   filter(!grepl("2012", IID))

# add 2012 labels 
year <-  ifelse(grepl("2012", raw_file$IID), "2012", "2014")
raw_file$FID <- paste0(raw_file$FID, "_", year)

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
  arrange(pop, chr, pos)

# write to file
write.table(pop_sums, file = "data/allele_freq_by_pop.txt", row.names = FALSE, quote = FALSE)

