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

plink_file <- list.files("data/snp_tables/", pattern = "outgroup_2014", full.names = TRUE)
output_folder <- "data/allele_frq"
#raw_file <- data.frame(fread(plink_file, header = TRUE, stringsAsFactors = FALSE))
raw_file <- read.table(plink_file[1], header = TRUE, stringsAsFactors = FALSE)

meta_df <- read.csv("metadata/mega_meta.csv")

# fix population names

#raw_file$FID <- gsub("NG-5241_.*", "DK", raw_file$FID)
#raw_file$FID <- gsub("SRX.*", "LC", raw_file$FID)
raw_file$FID <- gsub("[^A-Z]*", "", raw_file$FID)

# collapse CP pops
raw_file$FID[raw_file$FID == "CPN"|raw_file$FID == "CPSE"] <- "CP"
raw_file$FID[raw_file$FID == "BS"] <- "DK"
raw_file$IID <- gsub("BS", "DK", raw_file$IID)

# filter out low sample size pops
pop_counts <- meta_df %>% group_by(pop) %>% tally %>% filter(n > 4)
meta_df <- meta_df %>% 
  filter(pop %in% pop_counts$pop)

raw_file <- raw_file %>%
  filter(FID %in% pop_counts$pop)

# filter out 2012 data
raw_file <- raw_file %>% 
  filter(!grepl("2012|\\.2", IID))

meta_df$id <- gsub("BS", "DK", meta_df$id)

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

raw_long$FID <- gsub("SR_cmn", "GY_cmn", raw_long$FID)
raw_long$FID <- gsub("SR_wht", "GY_wht", raw_long$FID)

raw_long$FID <- gsub("SR_cmn|SF_cmn", "GY_cmn", raw_long$FID)
raw_long$FID <- gsub("SR_wht|SF_wht", "GY_wht", raw_long$FID)

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

# generate the allele ffrequencies
# arbitrarily, allele 0 

# merge target populations

pop_frq <- pop_sums %>%
  group_by(pop, chr, pos) %>%
  summarise(frq_a = sum(n[allele == 0]*2, n[allele == 1], na.rm = TRUE), frq_b = sum(n[allele == 2]*2, n[allele == 1], na.rm = TRUE)) %>%
  ungroup %>%
  mutate(frq_tot = frq_a + frq_b) %>%
  mutate(prop_a = frq_a / frq_tot, prop_b = frq_b / frq_tot)

pop_frq_b <- pop_frq %>%
  select(-frq_tot, -prop_a, -frq_a, -frq_b)

# spread the data and remove chr pos
pop_frq_b  <- spread(pop_frq_b, key = pop, value = prop_b) %>% 
  ungroup %>%
  arrange(chr, pos)

# select target populations
pop_frq_b <- pop_frq_b %>%
  select(chr, pos, SK_cbr, GY_cmn, GY_wht, DK_dk) %>% 
  filter(!is.na(DK_dk)) %>%
  filter(!is.na(SK_cbr)) %>%
  filter(!is.na(GY_cmn)) %>%
  filter(!is.na(GY_wht)) 

ab_df <- pop_frq_b %>%
  filter(DK_dk == 0) %>%
  filter(SK_cbr != 0) %>%
  filter(GY_cmn != 0) %>%
  filter(GY_wht != 0)


# frequency-based abba (from Martin et al. 2015 MBE)
# (1-p1)*p2*p3*(1-p4)
# if using p4 == 0 to find derived alleles, that term just goes to 1

ab_df$abba <- with(ab_df, (1-SK_cbr)*GY_cmn*GY_wht*(1-DK_dk))
ab_df$baba <- with(ab_df, SK_cbr*(1-GY_cmn)*GY_wht*(1-DK_dk))
ab_df$d_num <- with(ab_df, (abba - baba))
ab_df$d_den <- with(ab_df, (abba + baba))

d_stat <- sum(ab_df$d_num, na.rm = TRUE) / sum(ab_df$d_den, na.rm = TRUE)

ab_df %>%
  select(chr, pos, SK_cbr, GY_cmn, GY_wht, DK_dk) %>%
  gather(key = pop, value = frq, -chr, -pos) %>%
  mutate(pop_num = as.numeric(as.factor(pop))) %>%
  ggplot(aes(x = as.numeric(pos), y = frq + pop_num, color = pop))+
  geom_line()+
  facet_wrap(~chr)


# frequency-based baba (from Martin et al. 2015 MBE)
# p1*(1-p2)*p3*(1-p4)




pop_frq_a

#


write.table(pop_count[,1:2], "metadata/baypass_sites_2012.txt", row.names = FALSE, quote = FALSE)

pop_count <- pop_count %>%
  select(-chr, -pos)

# for some reason select isn't working, so this is a hack to fix that
#pop_count <- pop_count[,-1]

# filter for representation
#pop_count <- pop_count %>%
#  filter(!(DK_cmn %in% c("0,0", "1,1")))

#missing_matrix <- lapply(pop_count, function(x) as.numeric(x == "0,0" | is.na(x) | is.null (x))) %>% data.frame
#pop_count$missing <- rowSums(data.matrix(missing_matrix))

pop_count <- pop_count %>%
  #filter(missing < 2) %>%
  #select(-MH_cbr) %>%
  select(-SR_NA) %>%
  #select(-RT_cmn)

# convert to baypass format

pop_count_wht_cmn <- pop_count %>%
  select(matches("wht|cmn"))

pop_count_wht_cbr <-  pop_count %>%
  select(matches("wht|cbr"))

pop_count_cmn_cbr <-  pop_count %>%
  select(matches("cmn|cbr"))

split_bind <- function(x){
  strsplit(x %>% sapply(as.character), split = ',')  %>%
  do.call(rbind, .) %>% 
    data.frame
}

baypass_wht_cmn<- pop_count_wht_cmn %>% lapply(split_bind) %>% data.frame
baypass_wht_cbr <- pop_count_wht_cbr %>% lapply(split_bind) %>% data.frame
baypass_cmn_cbr <- pop_count_cmn_cbr %>% lapply(split_bind) %>% data.frame

#pop_count <- pop_count[,-16]
# out file names

base_name <- plink_file %>% 
  strsplit(split = "/") %>% 
  unlist %>% .[3]

out_file_baypass_wht_cmn <- base_name %>%
  gsub(".gz", ".baypass.wht_cmn.2012.geno", .) %>%
  paste0(output_folder, "/", .)

out_file_baypass_wht_cbr <- base_name %>%
  gsub(".gz", ".baypass.wht_cbr.geno", .) %>%
  paste0(output_folder, "/", .)

out_file_baypass_cmn_cbr <- base_name %>%
  gsub(".gz", ".baypass.cmn_cbr.geno", .) %>%
  paste0(output_folder, "/", .)
  

# write to file
write.table(baypass_wht_cmn, file = out_file_baypass_wht_cmn, row.names = FALSE, quote = FALSE, col.names = FALSE)
write.table(baypass_wht_cbr, file = out_file_baypass_wht_cbr, row.names = FALSE, quote = FALSE, col.names = FALSE)
write.table(baypass_cmn_cbr, file = out_file_baypass_cmn_cbr, row.names = FALSE, quote = FALSE, col.names = FALSE)

# write baypass "ecotype" files

wht_cmn_efile <- names(pop_count_wht_cmn) %>%
  gsub("[A-Z]{2}_", "", .) %>%
  gsub("wht", "1", .) %>%
  gsub("cmn", "-1", .)

wht_cbr_efile <- names(pop_count_wht_cbr) %>%
  gsub("[A-Z]{2}_", "", .) %>%
  gsub("wht", "1", .) %>%
  gsub("cbr", "-1", .) 

cmn_cbr_efile <- names(pop_count_cmn_cbr) %>%
  gsub("[A-Z]{2}_", "", .) %>%
  gsub("cmn", "1", .) %>%
  gsub("cbr", "-1", .) 


write.table(t(wht_cmn_efile), file = gsub("geno", "efile", out_file_baypass_wht_cmn), row.names = FALSE, quote = FALSE, col.names = FALSE)
write.table(t(wht_cbr_efile), file = gsub("geno", "efile", out_file_baypass_wht_cbr), row.names = FALSE, quote = FALSE, col.names = FALSE)
write.table(t(cmn_cbr_efile), file = gsub("geno", "efile", out_file_baypass_cmn_cbr), row.names = FALSE, quote = FALSE, col.names = FALSE)




