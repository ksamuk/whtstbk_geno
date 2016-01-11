#### fix the malformed FASTA headers from PGDSpider VCF -> FASTA

library(dplyr)
library(stringr)
library(ape)

# find the fasta file
fasta_file <- list.files("data/other_formats", pattern = "fasta.fa$", full.names = TRUE)

# read in as dnaBIN (see ape docs)
fasta <- read.dna(fasta_file, format = "fasta", skip = 0, as.character = TRUE, as.matrix = FALSE)

# grab the headers
fasta_header <- names(fasta)

# WIELD THE MIGHTY HAMMER 
fixed_headers <- fasta_header %>% 
  as.character %>% 
  gsub("\\|.*\\| *|locus:[12]*|SNP_[0-9]*|,", "", .) %>% 
  gsub("whtstbk_gbs_2012_brds_", "2012_", .) %>%
  gsub("\\.", "_", .) %>%
  str_trim

# apply fix  
names(fasta) <- fixed_headers

replace_n <- function(x){ 
  x <- toupper(x)
  x[x=="N"]<-"-"
  x
  }

fasta_out <- lapply(fasta, replace_n)

# drop 1/2 of sequences (read 1 and 2, apply IUPAC genotype names later)

file_name <- fasta_file %>% gsub(".fa$", ".fixed.fa", .)
write.dna(fasta_out, file_name, format = "fasta", nbcol = -6, colsep = "")
#write.nexus.data(fasta, file = file_name)
