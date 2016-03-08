# run rehh on the output of a shapeit file (.haps)

################################################################################
# Libraries
################################################################################

rm(list=ls())
library("rehh")
library("dplyr")
library("ggplot2")

list.files("functions", full.names = TRUE) %>% sapply(.,source, verbose = FALSE, echo = FALSE) %>% invisible

################################################################################
# Input file locations
################################################################################

# input file locations

# the .out files from fastphase
out.files <- list.files("data/phased_fastphase", pattern = ".out", full.names = TRUE)
map.files  <- list.files("data/map_files", pattern = ".map", full.names = TRUE)

# meta data
meta_df <- read.csv("metadata/mega_meta.csv")

################################################################################
# Build .map files
################################################################################

chr_set <- as.roman(rep(1:21)) %>% as.character

map_master <- read.table("metadata/ref_alt_2014_bial_maf.txt", header = TRUE, stringsAsFactors = FALSE)

build_map_file <- function(chr, folder = "data/phased_fastphase/input_maps"){
  
  hap_file <- list.files(folder, pattern = paste0("^","chr",chr,".map"), full.names = TRUE)
  map_file <- read.table(hap_file, header = FALSE)[,c(2,4)]
  
  map_chr <- map_file[,1] %>% gsub(":[0-9]*", "", .)
  map_file <- data.frame(chr = map_chr, pos = map_file[,2])
  map_file <- left_join(map_file, map_master, by = c("chr", "pos"))
  
  snp_name <- paste0(map_file$chr, "_", map_file$pos)
  map_file <- cbind(snp_name, map_file)
  
  names(map_file)[2:5] <- c("chr", "pos", "anc", "der")
  
  map_file$anc <- 0
  map_file$der <- 1
  
  map_file$chr <- map_file$chr %>% gsub("chr", "", .) %>% as.roman %>% as.numeric
  
  file_name <- paste0("data/map_files/", "chr", chr, ".map")
  
  write.table(map_file, file = file_name, row.names = FALSE, quote = FALSE, col.names = FALSE)
  
}

lapply(chr_set, build_map_file)


################################################################################
# Parse fastphase output
################################################################################


fp_pop <- read.table("metadata/pop_file_fastphase.txt", header = FALSE)
names(fp_pop) <- c("id", "pop", "cluster")
fp_pop$number <- as.numeric(fp_pop[,2])

merge_fastphase_output <- function(chr, species, fp_pop = fp_pop){
  
  out_file <- grep(paste0("chr", chr,"_"), out.files, value = TRUE)
  map_file <- grep(paste0("chr", chr,".map"), map.files, value = TRUE)
  
  pop_nums <- fp_pop %>%
    filter(cluster == species) %>%
    select(number) %>%
    unique %>% unlist %>% as.numeric
  
  # combine all the individual haplotypes into a single .hh file
  
  haplohh_cluster <- data2haplohh(out_file, map_file, popsel = pop_nums[1], recode.allele = TRUE)
  
  for (i in 2:length(pop_nums)){
    haplohh_tmp <- data2haplohh(out_file, map_file, popsel = pop_nums[i], recode.allele = TRUE)
    
    slot(haplohh_cluster, "haplo") <- rbind(slot(haplohh_cluster, "haplo"), slot(haplohh_tmp, "haplo"))  
    slot(haplohh_cluster, "nhap") <- sum(slot(haplohh_cluster, "nhap"), slot(haplohh_tmp, "nhap")) 
    
  }
  
  return(haplohh_cluster)
}

create_rehh_ehh_file <- function(chr, meta_df, haps_dir = "data/phased_fastphase"){
  
  cmn_hh <- merge_fastphase_output(chr, "cmn", fp_pop)
  wht_hh <- merge_fastphase_output(chr, "wht", fp_pop)
  cbr_hh <- merge_fastphase_output(chr, "cbr", fp_pop)
  
  # read in an perform hh scan
  cat ("Performing hh scans...")
  
  hh_scans <- list()
  ihs_scans <- list()
  hh_scans_rsb <- list()
  
  for (j in 1:3){
    
    i <- c("wht", "cmn", "cbr")[j]
    hh_scans[[j]] <- scan_hh(get(paste0(i, "_hh")))
    ihs_scans[[j]] <- ihh2ihs(hh_scans[[j]])$res.ihs
    
    ihs_scans[[j]] <- ihs_scans[[j]]%>% data.frame %>%
      dplyr::select(CHR, POSITION, iHS)
    names(ihs_scans[[j]]) <- c("chr", "pos", paste0("ihs_", i))
    
    hh_scans_rsb[[j]] <- hh_scans[[j]]
    hh_scans[[j]] <- hh_scans[[j]]%>% data.frame %>%
      dplyr::select(CHR, POSITION, IES)
    
    names(hh_scans[[j]]) <- c("chr", "pos", paste0("ies_", i))
    #save(ehh.scan, file = paste0("data/rehh/ehh_scans/chr",chr,".ehh.",i,".Rdata")
    
  }
  
  ies_df <- left_join(hh_scans[[1]], hh_scans[[2]]) %>% left_join(hh_scans[[3]])
  ihs_df <- left_join(ihs_scans[[1]], ihs_scans[[2]]) %>% left_join(ihs_scans[[3]])
  
  write.table(ies_df, paste0("data/rehh/ehh_scans/chr",chr,".txt"), row.names = FALSE, quote = FALSE)
  write.table(ihs_df, paste0("data/rehh/ihs_scans/chr",chr,".txt"), row.names = FALSE, quote = FALSE)
  
  hh_scans <- hh_scans_rsb
  
  cat ("Performing Rsb scans...")
  
  wht_cmn_rsb <- ies2rsb(hh_scans[[1]], hh_scans[[2]], popname1 = "wht", popname2 = "cmn", method="bilateral") %>% data.frame
  names(wht_cmn_rsb) <- c("chr", "pos", "rsb_wht_cmn", "rsb_pval_wht_cmn")
  
  wht_cbr_rsb <- ies2rsb(hh_scans[[1]], hh_scans[[3]], popname1 = "wht", popname2 = "cbr", method="bilateral") %>% data.frame
  names(wht_cbr_rsb) <- c("chr", "pos", "rsb_wht_cbr", "rsb_pval_wht_cbr")
  
  cbr_cmn_rsb <- ies2rsb(hh_scans[[3]], hh_scans[[2]], popname1 = "cbr", popname2 = "cmn", method="bilateral") %>% data.frame
  names(cbr_cmn_rsb) <- c("chr", "pos", "rsb_cmn_cbr", "rsb_pval_cmn_cbr")
  
  rsb_df <- left_join(wht_cmn_rsb, wht_cbr_rsb, by = c("chr", "pos")) %>% left_join(cbr_cmn_rsb, by = c("chr", "pos"))
  head(rsb_df)
  
  write.table(rsb_df, paste0("data/rehh/rsb_scans/chr",chr,".txt"), row.names = FALSE, quote = FALSE)
  

}


################################################################################
# Release the hounds!
################################################################################

chr_set <- as.roman(rep(1:21)) %>% as.character

# calc ehh and rsb for all the things
lapply(chr_set, create_rehh_ehh_file, meta_df = meta_df)

# create master ehh and rsb files

ehh.files <- list.files("data/rehh/ehh_scans", full.names = TRUE)
ehh_df <- lapply(ehh.files, read.table, header = TRUE)
ehh_df <- bind_rows(ehh_df)
write.table(ehh_df, "data/stats/whtstbk_ies_df.txt", quote = FALSE, row.names = FALSE)

rsb.files <- list.files("data/rehh/rsb_scans", full.names = TRUE)
rsb_df <- lapply(rsb.files, read.table, header = TRUE)
rsb_df <- bind_rows(rsb_df)
write.table(rsb_df , "data/stats/whtstbk_rsb_df.txt", quote = FALSE, row.names = FALSE)

ihs.files <- list.files("data/rehh/ihs_scans", full.names = TRUE)
ihs_df <- lapply(ihs.files, read.table, header = TRUE)
ihs_df <- bind_rows(ihs_df)
write.table(ihs_df , "data/stats/whtstbk_ihs_df.txt", quote = FALSE, row.names = FALSE)

