# format a .raw plink file to treemix AND BayPass format
# assumes plink is malformed from vcf conversion
# e.g. the FID and IID columns are identical

# Run like: Rscript format_raw_to_treemix_baypass.R /path/to/zipped/.raw/file /output/folder

#library("data.table")

#library("readr")
#library("data.table")
#library("stringr")
library("tidyr")
library("dplyr")

list.files("functions", full.names = TRUE) %>% sapply(.,source, verbose = FALSE, echo = FALSE) %>% invisible
select <- dplyr::select

args <- commandArgs(trailingOnly = TRUE)
#plink_file <- args[1]
plink_file <- list.files("data/snp_tables/", pattern = "outgroup_2014", full.names = TRUE)
#output_folder <- args[2]
output_folder <- "data/treemix"

#raw_file <- data.frame(fread(plink_file, header = TRUE, stringsAsFactors = FALSE))
raw_file <- read.table(plink_file, header = TRUE, stringsAsFactors = FALSE)

meta_df <- read.csv("metadata/mega_meta.csv")

# fix population names

raw_file$FID <- gsub("SRX.*", "LC", raw_file$FID)
raw_file$FID <- gsub("BS*", "DK", raw_file$FID)
raw_file$FID <- gsub("[^A-Z]*", "", raw_file$FID)

# filter out low sample size pops
# * this might be redundant in new workflow *

# collapse CP pops
raw_file$FID[raw_file$FID == "CPN"|raw_file$FID == "CPSE"] <- "CP"

pop_counts <- meta_df %>% group_by(pop) %>% tally %>% filter(n > 4)
meta_df <- meta_df %>% 
  filter(pop %in% pop_counts$pop)

raw_file <- raw_file %>%
  filter(FID %in% pop_counts$pop)

# filter out 2012 data
raw_file <- raw_file %>% 
  filter(!grepl("2012|\\.2", IID))

# add in species labels
pop_suffixes <- left_join(data.frame(id = raw_file$IID), meta_df)
raw_file$FID <- paste0(pop_suffixes$pop, "_", pop_suffixes$cluster)

# remove all but genotype data
raw_file <- raw_file %>%
  .[, -c(2:6)]

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
  ungroup %>%
  arrange(chr, pos) 

write.table(pop_count[,1:2], "metadata/treemix_sites.txt", row.names = FALSE, quote = FALSE)

pop_count <- pop_count %>%
  select(-chr, -pos)

# for some reason select isn't working, so this is a hack to fix that
#pop_count <- pop_count[,-1]

# filter for representation
#pop_count <- pop_count %>%
#  filter(!(DK_cmn %in% c("0,0", "1,1")))

#missing_matrix <- lapply(pop_count, function(x) as.numeric(x == "0,0" | is.na(x) | is.null (x))) %>% data.frame
#pop_count$missing <- rowSums(data.matrix(missing_matrix))

# RT_cmn is pushing it

pop_count <- pop_count %>%
  #filter(missing < 2) %>%
  #select(-MH_cbr) %>%
  select(-SR_NA) %>%
  select(-LN_NA) %>%
  select(-RT_cmn)

#pop_count <- pop_count[,-16]
# out file names

base_name <- plink_file %>% 
  strsplit(split = "/") %>% 
  unlist %>% .[3]

out_file_treemix <- base_name %>%
  gsub(".gz", ".treemix.gz", .) %>%
  paste0(output_folder, "/", .)

# write to file
write.table(pop_count, file = gzfile(out_file_treemix), row.names = FALSE, quote = FALSE)





