# format a .raw plink file to treemix AND BayPass format
# assumes plink is malformed from vcf conversion
# e.g. the FID and IID columns are identical

# Run like: Rscript format_raw_to_baypass.R /path/to/zipped/.raw/file /output/folder

library("tidyr")
library("dplyr")


args <- commandArgs(trailingOnly = TRUE)
plink_file <- args[1]
output_folder <- args[2]

list.files("functions", full.names = TRUE) %>% sapply(.,source, verbose = FALSE, echo = FALSE) %>% invisible

select <- dplyr::select

plink_file <- list.files("data/snp_tables/", pattern = "whtstbk_bial_nomaf_nosex.raw", full.names = TRUE)
output_folder <- "data/dadi"
#raw_file <- data.frame(fread(plink_file, header = TRUE, stringsAsFactors = FALSE))
raw_file <- data.frame(data.table::fread(plink_file[1]))
meta_df <- read.csv("metadata/mega_meta.csv")

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


snp_chr_to_dadi <- function(raw_file){
  
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
  
  # collapse by *species (cluster)* label (filter out CBR)
  raw_file <- raw_file[!grepl("cbr", raw_file$FID),]
  raw_file$FID <- ifelse(grepl("wht", raw_file$FID), "wht", "cmn")
  
  
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
    arrange(chr, pos)
  
  pop_counts <- pop_sums %>%
    group_by(pop, chr, pos) %>%
    summarise(frq_a = sum(n[allele == 0]*2, n[allele == 1], na.rm = TRUE), frq_b = sum(n[allele == 2]*2, n[allele == 1], na.rm = TRUE)) %>%
    ungroup
  
  # filter invariant sites
  # these can be retained due to their presence in other populations in the snp table that were filtered out earlier
  pop_counts <- pop_counts %>%
    arrange(chr, pos) %>%
    group_by(chr, pos) %>%
    mutate(variant = sum(frq_a) > 0 && sum(frq_b) > 0) %>%
    ungroup %>%
    filter(variant == TRUE) %>%
    select(-variant)
  
  # spread the data
  # one data frame per allele
  pop_counts_0  <- pop_counts %>%
    select(-frq_b) %>%
    spread(key = pop, value = frq_a) %>% 
    arrange(chr, pos)
  
  pop_counts_0$pos <- as.numeric(as.character(pop_counts_0$pos))
  names(pop_counts_0)[3:4] <- paste0(names(pop_counts_0)[3:4], "1")
  
  pop_counts_1  <- pop_counts %>%
    select(-frq_a) %>%
    spread(key = pop, value = frq_b) %>% 
    arrange(chr, pos)
  
  pop_counts_1$pos <- as.numeric(as.character(pop_counts_1$pos))
  names(pop_counts_1)[3:4] <- paste0(names(pop_counts_1)[3:4], "2")
  
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
    dplyr::select(ref, alt, Allele1, cmn1, wht1, Allele2, cmn2, wht2, chr, pos)
  
  
  dadi_df
  
}

# remove datasets with no sites (?!?)

dadi_list <- lapply(snp_chr, snp_chr_to_dadi)

dadi_df <- bind_rows(dadi_list) %>%
  arrange(chr, pos)

names(dadi_df)[4:5] <- c("cmn", "wht")
names(dadi_df)[7:8] <- c("cmn", "wht")

write.table(dadi_df, "data/dadi/whtstbk_2014_nomaf_bial_nosex.txt", quote = FALSE, row.names = FALSE)

dadi_df <- read.table("data/dadi/whtstbk_2014_nomaf_bial_nosex.txt", header = TRUE)

dadi_df <- dadi_df %>% 
  mutate(sum_counts = cmn + wht + cmn.1 + wht.1) %>%
  filter(sum_counts > 350) %>%
  dplyr::select(-sum_counts)

dadi_df <- dadi_df %>% 
  dplyr::select(alt, ref, Allele1, cmn, wht, Allele2, cmn.1, wht.1, chr, pos)

names(dadi_df)[4:5] <- c("cmn", "wht")
names(dadi_df)[7:8] <- c("cmn", "wht")

write.table(dadi_df, "data/dadi/whtstbk_2014_nomaf_bial_nosex_filtered.txt", quote = FALSE, row.names = FALSE)
  
dadi_df <- dadi_df %>% 
  filter(cmn != 1, wht != 1, cmn.1 != 1, wht.1 !=1)

write.table(dadi_df, "data/dadi/whtstbk_2014_no_singletons.txt", quote = FALSE, row.names = FALSE)

