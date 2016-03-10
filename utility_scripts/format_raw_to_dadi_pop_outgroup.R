# format a .raw plink file to treemix AND BayPass format
# assumes plink is malformed from vcf conversion
# e.g. the FID and IID columns are identical

# Run like: Rscript format_raw_to_dadi.R /path/to/zipped/.raw/file /output/folder

################################################################################
# Libraries
################################################################################

library("tidyr")
library("dplyr")

list.files("functions", full.names = TRUE) %>% sapply(.,source, verbose = FALSE, echo = FALSE) %>% invisible
select <- dplyr::select

################################################################################
# Initials / Inputs
################################################################################

# parse args (if any)
args <- commandArgs(trailingOnly = TRUE)
plink_file <- args[1]
output_folder <- args[2]

# the main inputs
plink_file <- list.files("data/snp_tables/", pattern = "whtstbk_bial_nomaf_nosex_outgroup_pruned.raw", full.names = TRUE)
#plink_file <- list.files("data/snp_tables/", pattern = "whtstbk_2014_pruned.raw", full.names = TRUE)

output_folder <- "data/dadi"

# read in raw file, meta data and sites to exclude (e.g. selected sites)
raw_file <- data.frame(data.table::fread(plink_file[1]))
meta_df <- read.csv("metadata/mega_meta.csv")
non_outliers_df <- read.table("metadata/non_outlier_sites.txt", h = T)
non_outliers_df <- non_outliers_df %>%
  mutate(chr_pos = paste0(chr, "_", pos))

meta_df$cluster[meta_df$pop == "DK"] <- "cmn"

meta_df$pop <- as.character(meta_df$pop)
meta_df$pop[meta_df$pop == "DK"] <- "BS"

################################################################################
# Processing
################################################################################

# split the raw file by chromosome
col_names <- names(raw_file)
chr_names <- col_names[-c(1:6)] %>% gsub("\\..*", "", .) %>% unique

extract_chr_columns_snp_table <- function(chr_name, col_names, raw_file){
  
  non_geno <- raw_file[,1:6]
  chr_cols <- grepl(paste0(chr_name, "\\."), col_names)
  raw_sub <- raw_file[,chr_cols]
  cbind(non_geno, raw_sub)

}

snp_chr <- lapply(chr_names, extract_chr_columns_snp_table, raw_file = raw_file, col_names = col_names)
snp_chr <- snp_chr[lapply(snp_chr, length) > 6]

rm(raw_file)

snp_chr_to_dadi <- function(raw_file, non_outliers_df = non_outliers_df){
  
  cat(".")
  
  # fix population names
  raw_file$FID <- gsub("[^A-Z]*", "", raw_file$FID)
  
  # collapse CP pops
  raw_file$FID[raw_file$FID == "CPN"|raw_file$FID == "CPSE"] <- "CP"
  
  # filter out 2012 data
  raw_file <- raw_file %>% 
    filter(!grepl("2012|\\.2", IID))
  
  # filter out low sample size pops
  pop_counts <- meta_df %>% group_by(pop) %>% tally %>% filter(n > 4)
  meta_df <- meta_df %>% 
    filter(pop %in% pop_counts$pop)
  
  raw_file <- raw_file %>%
    filter(FID %in% pop_counts$pop)
  
  # add in species labels
  pop_suffixes <- left_join(data.frame(id = raw_file$IID), meta_df, by = "id")
  raw_file$FID <- paste0(pop_suffixes$pop, "_", pop_suffixes$cluster)
  
  # remove all but genotype data
  raw_file <- raw_file %>%
    .[, -c(2:6)]
  
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
    arrange(chr, pos) %>%
    filter(!grepl("NA", pop))
  
  pop_counts <- pop_sums %>%
    group_by(pop, chr, pos) %>%
    summarise(frq_a = sum(n[allele == 0]*2, n[allele == 1], na.rm = TRUE), frq_b = sum(n[allele == 2]*2, n[allele == 1], na.rm = TRUE)) %>%
    ungroup
  
  # filter sites to exclude outliers (non_outliers_df)
  
  pop_counts <- pop_counts %>%
    mutate(chr_pos = paste0(chr, "_", pos)) %>%
    filter(chr_pos %in% non_outliers_df$chr_pos) %>%
    select(-chr_pos)
  
  # filter invariant sites
  # these can be retained due to their presence in other populations in the snp table that were filtered out earlier
  #pop_counts <- pop_counts %>%
  #  arrange(chr, pos) %>%
  #  group_by(chr, pos) %>%
  #  mutate(variant = sum(frq_a) > 0 | sum(frq_b) > 0) %>%
  #  ungroup %>%
  #  filter(variant == TRUE) %>%
  #  select(-variant)
  
  # spread the data
  # one data frame per allele
  pop_counts_0  <- pop_counts %>%
    select(-frq_b) %>%
    spread(key = pop, value = frq_a) %>% 
    arrange(chr, pos)
  
  pop_counts_0$pos <- as.numeric(as.character(pop_counts_0$pos))
  names(pop_counts_0)[-c(1:2)] <- paste0(names(pop_counts_0)[-c(1:2)], "1")
  
  pop_counts_1  <- pop_counts %>%
    select(-frq_a) %>%
    spread(key = pop, value = frq_b) %>% 
    arrange(chr, pos)
  
  pop_counts_1$pos <- as.numeric(as.character(pop_counts_1$pos))
  names(pop_counts_1)[-c(1:2)] <- paste0(names(pop_counts_1)[-c(1:2)], "2")
  
  # join in ref info
  
  ref_alt <- read.table("metadata/ref_alt_nova_scotia_master.txt", header = TRUE, stringsAsFactors = FALSE)
  ref_alt$chr <- ref_alt$chr %>% gsub("chr", "", .) %>% gsub("Un", "XXII", .) %>% as.roman %>% as.numeric
  
  dadi_df <- left_join(pop_counts_0, pop_counts_1, by = c("chr", "pos"))
  dadi_df <- left_join(dadi_df, ref_alt, by = c("chr", "pos"))
  
  dadi_df$Allele1 <- dadi_df$ref
  dadi_df$Allele2 <- dadi_df$alt
  
  dadi_df$ref <- paste0("-", dadi_df$ref, "-")
  dadi_df$alt <- paste0("-", dadi_df$alt, "-")
  
  dadi_df <- dadi_df %>%
    dplyr::select(alt, ref, Allele1, matches("1$", ignore.case = TRUE), Allele2, matches("2$", ignore.case = TRUE), chr, pos)
  
  
  dadi_df
  
}

# apply funciton to chr list
dadi_list <- lapply(snp_chr, snp_chr_to_dadi, non_outliers_df = non_outliers_df)

# bind rows
dadi_df <- bind_rows(dadi_list) %>%
  arrange(chr, pos)

# write master to file
write.table(dadi_df, "data/dadi/whtstbk_dk_dadi_pops_master.txt", quote = FALSE, row.names = FALSE)

#write.table(dadi_df, "data/dadi/whtstbk_2014_dadi_pops_master.txt", quote = FALSE, row.names = FALSE)

# can manipulate master staritng here
dadi_df <- read.table("data/dadi/whtstbk_2014_dadi_pops_master.txt", header = TRUE)

# coverage filter (drop low cov loci)
counts_df <- dadi_df %>% 
  dplyr::select(-alt, -ref, -Allele1, -Allele2, -chr, -pos) %>%
  mutate(count = rowSums(.))

dadi_df <- data.frame(dadi_df, count = counts_df$count)

dadi_df <- dadi_df %>%
  filter(count < 540) %>%
  dplyr::select(-count)

dadi_df <- dadi_df%>% 
  filter(!is.na(Allele1))

# make an "analysis ready" file (with forced duplicate column names)
len_names <- length(names(dadi_df))

# god forgive me for my failing regex here
# its beacuse of "allele1" and "allele2"
names(dadi_df)[-c(1:2, (len_names-1):len_names)] <- names(dadi_df)[-c(1:2, (len_names-1):len_names)] %>% 
  gsub("cmn1", "cmn", .) %>% 
  gsub("cbr1", "cbr", .) %>% 
  gsub("wht1", "wht", .) %>% 
  gsub("cmn2", "cmn", .) %>% 
  gsub("cbr2", "cbr", .) %>% 
  gsub("wht2", "wht", .) 

write.table(dadi_df, "data/dadi/whtstbk_2014_dadi_pops_outgroup_analysis.txt", quote = FALSE, row.names = FALSE)

# merge some populations of choice
# e.g. because they are really close geographically

wht_pops_to_merge <- c("SR_wht", "SF_wht", "MH_wht")
cmn_pops_to_merge <- c("SR_cmn", "SF_cmn", "MH_cmn")

dadi_sub <- dadi_df %>%
  select(alt, ref, Allele1, Allele2, chr, pos)

dadi_wht <- dadi_df %>%
  select(contains("SR_wht"), contains("SF_wht"), contains("MH_wht"))

dadi_cmn <- dadi_df %>%
  select(contains("SR_cmn"), contains("SF_cmn"), contains("MH_cmn"))

dadi_wht <- dadi_wht %>%
  mutate(wht1 = SR_wht1 + SF_wht1 + MH_wht1) %>%
  mutate(wht2 = SR_wht2 + SF_wht2 + MH_wht2) %>%
  select(wht1, wht2)

dadi_cmn <- dadi_cmn %>%
  mutate(cmn1 = SR_cmn1 + SF_cmn1 + MH_cmn1) %>%
  mutate(cmn2 = SR_cmn2 + SF_cmn2 + MH_cmn2) %>%
  select(cmn1, cmn2)

dadi_sub <- data.frame(dadi_sub[,1:3], dadi_wht$wht1, dadi_cmn$cmn1, dadi_sub$Allele2, dadi_wht$wht2, dadi_cmn$cmn2, dadi_sub$chr, dadi_sub$pos)
dadi_sub <- dadi_sub%>% 
  filter(!is.na(Allele1))

names(dadi_sub)[4:10] <- c("wht", "cmn", "Allele2", "wht", "cmn", "chr", "pos")

write.table(dadi_sub , "data/dadi/whtstbk_2014_dadi_pops_collapsed.txt", quote = FALSE, row.names = FALSE)

# drop singletons (seq errors) 
dadi_df <- dadi_df %>% 
  filter(cmn != 1, wht != 1, cmn.1 != 1, wht.1 !=1)

write.table(dadi_df, "data/dadi/whtstbk_2014_no_singletons.txt", quote = FALSE, row.names = FALSE)

