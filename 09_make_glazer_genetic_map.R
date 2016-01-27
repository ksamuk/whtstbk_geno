#### Subsets, filters and reformats a VCF 
#### end point is a 'geno_df' 
#### example run: Rscript 01_process_vcf_to_geno_df.R data/vcf/whtstbk_master.vcf.bgz multi 2

################################################################################
# Libraries
################################################################################

library("dplyr")
library("readxl")
library("stringr")

################################################################################
# Data input and formatting
################################################################################

# read in glazer et al orig vs. revised mappings

pos_dat <- read_excel("metadata/orig_revised_positions_glazer.xlsx")

ftc_map <- read_excel("metadata/ftc_genetic_map_glazer.xlsx")
bepa_map <- read_excel("metadata/bepa_genetic_map_glazer.xlsx")
bepa_map <- bepa_map[,1:3]

#format and join physical and map position data
new_pos_map <- pos_dat %>%
  select(newChr, marker, newStart, newEnd)

names(new_pos_map) <- c("chr", "marker", "start", "end")
names(ftc_map) <- c("chr", "marker", "gen_pos_ftc")
names(bepa_map) <- c("chr", "marker", "gen_pos_bepa")

gen_pos <- left_join(new_pos_map, ftc_map)
gen_pos <- left_join(gen_pos, bepa_map)

################################################################################
# Calculate recombination rates (cM/MB)
################################################################################

gen_pos <- gen_pos %>%
  rowwise %>%
  mutate(mid = floor(mean(c(start, end)))) %>%
  mutate(start_new = ifelse(start < end, start, end), end_new = ifelse(start < end, end, start)) %>%
  mutate(start = start_new, end = end_new) %>%
  select(-start_new, -end_new) %>%
  arrange(chr, start)

recomb_rates <- gen_pos %>%
  select(chr, start, end, mid, gen_pos_ftc, gen_pos_bepa)

recomb_rates <- recomb_rates %>%
  select(-start, -end) %>%
  group_by(chr) %>%
  mutate(pos1 = mid, pos2 = lead(mid)) %>%
  mutate(phys_dist = lead(mid) - mid) %>%
  mutate(ftc_dist = lead(gen_pos_ftc) - gen_pos_ftc)%>%
  mutate(bepa_dist = lead(gen_pos_bepa) - gen_pos_bepa)

View(recomb_rates)


# 





