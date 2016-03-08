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
  mutate(ftc_dist = abs(lead(gen_pos_ftc) - gen_pos_ftc))%>%
  mutate(bepa_dist = abs(lead(gen_pos_bepa) - gen_pos_bepa)) %>%
	ungroup

recomb_rates <- recomb_rates %>%
	filter(!is.na(pos2)) %>%
	filter(!is.na(gen_pos_ftc) | !is.na(gen_pos_bepa)) %>%
	rowwise %>%
	mutate(avg_dist = mean(c(ftc_dist, bepa_dist), na.rm = TRUE)) %>%
	mutate(recomb_rate = avg_dist / phys_dist) %>%
	rowwise %>%
	mutate(mid = floor(mean(c(pos1, pos2)))) %>%
	select(chr, mid, recomb_rate, avg_dist) %>%
	filter(!is.nan(recomb_rate))

# just FTC maps
# recomb_rates_ftc <- recomb_rates %>%
# 	select(-contains("bepa")) %>%
# 	filter(!is.na(pos2)) %>%
# 	filter(!is.na(gen_pos_ftc)) %>%
# 	mutate(recomb_rate = ftc_dist / phys_dist) %>%
# 	select(chr, mid, recomb_rate, gen_pos_ftc) %>%
# 	filter(!is.na(recomb_rate))


	
write.table(recomb_rates, file = "metadata/shapeit_maps/glazer_recomb.txt", quote = FALSE, row.names = FALSE)



