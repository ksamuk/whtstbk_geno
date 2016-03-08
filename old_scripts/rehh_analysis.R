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

haps.files <- list.files("data/phased", pattern = ".haps", full.names = TRUE)
sample.files <- list.files("data/phased", pattern = ".sample", full.names = TRUE)
meta_df <- read.csv("metadata/mega_meta.csv")


################################################################################
# Release the hounds!
################################################################################

create_rehh_ehh_file <- function(chr, meta_df, haps_dir = "data/phased"){
  
  #Map file contains SNPs information in five columns SNP names, chromosome, position,
  #ancestral and derived allele.
  cat(paste0(chr,"..."))
  cat("Formatting phased hap files...")
  
  #make the map file
  
  haps.file <- list.files("data/phased", pattern = paste0("chr",chr,"\\..*haps"), full.names = TRUE)
  haps.file <- read.table(haps.file, stringsAsFactors = FALSE)
  sample.file <- list.files("data/phased", pattern = paste0("chr",chr,"\\..*samp"), full.names = TRUE)
  sample.file <- read.table(sample.file, header = TRUE, stringsAsFactors = FALSE)
  sample.file <- sample.file[-1,]
  
  map.file <- data.frame(snp.name = 1:length(haps.file[,1]) , lg = chrom.to.num(haps.file[,1]), pos = haps.file[,3], anc = rep(1, length(haps.file[,1])) , der = rep(2, length(haps.file[,1])))
  
  # make the hap input file
  
  hap.input <- haps.file[,-(1:5)] %>% as.matrix
  row.names(hap.input) <- haps.file[,3]
  colnames(hap.input) <- rep(sample.file$ID_1, each = 2)
  hap.input <- hap.input + 1
  hap.input <- t(hap.input)
  
  #annotate rows properly (no dup rows allows)
  
  haplo.names <- rep(sample.file$ID_1, each = 2)
  first.haplos <- seq(from = 1, to = length(haplo.names), by = 2)
  second.haplos <- seq(from = 2, to = length(haplo.names), by = 2)
  
  haplo.names[first.haplos] <- paste0(haplo.names[first.haplos],"_1")
  haplo.names[second.haplos] <- paste0(haplo.names[second.haplos],"_2")
  
  # get the cluster assignments for each id (for splitting later)
  short_names <- data.frame(id = gsub("_[12]", "", haplo.names))
  short_names <- left_join(short_names, meta_df, by = "id") %>% 
    select(id, cluster)
  
  # the master hap.input file
  hap.input <- data.frame(cluster = short_names$cluster, haplo.names, hap.input)
  
  hap_wht <- hap.input %>% filter(cluster == "wht") %>% select(-cluster)
  hap_cmn <- hap.input %>% filter(cluster == "cmn") %>% select(-cluster)
  hap_cbr <- hap.input %>% filter(cluster == "cbr") %>% select(-cluster)
  
  # write to file (necessary for rehh for some reason)
  # the map file (same for all pops)
  write.table(map.file, file = paste0("data/rehh/map/chr", chr, ".map"), quote = FALSE, col.names = FALSE, row.names = FALSE)
  
  # the .hap files x 3
  write.table(hap_wht, file = paste0("data/rehh/wht_hap/chr", chr, ".wht.hap"), quote = FALSE, col.names = FALSE, row.names = FALSE)
  write.table(hap_cmn, file = paste0("data/rehh/cmn_hap/chr", chr, ".cmn.hap"), quote = FALSE, col.names = FALSE, row.names = FALSE)
  write.table(hap_cbr, file = paste0("data/rehh/cbr_hap/chr", chr, ".cbr.hap"), quote = FALSE, col.names = FALSE, row.names = FALSE)
  
  # read in an perform hh scan
  cat ("Performing hh scans...")
  
  hh_scans <- list()
  ihs_scans <- list()
  hh_scans_rsb <- list()
  
  for (j in 1:3){
    
    i <- c("wht", "cmn", "cbr")[j]
    hap.in <- data2haplohh(hap_file = paste0("data/rehh/",i,"_hap/chr", chr, ".",i,".hap"), map_file = paste0("data/rehh/map/chr", chr, ".map"))
    
    hh_scans[[j]] <- scan_hh(hap.in)
    ihs_scans[[j]] <- ihh2ihs(hh_scans[[j]])$res.ihs 
    hh_scans_rsb[[j]] <- hh_scans[[j]]
    hh_scans[[j]] <- hh_scans[[j]]%>% data.frame %>%
      dplyr::select(CHR, POSITION, IES)
    
    names(hh_scans[[j]]) <- c("chr", "pos", paste0("ies_", i))
    #save(ehh.scan, file = paste0("data/rehh/ehh_scans/chr",chr,".ehh.",i,".Rdata")
    
  }
  
  ies_df <- left_join(hh_scans[[1]], hh_scans[[2]]) %>% left_join(hh_scans[[3]])
  
  write.table(ies_df, paste0("data/rehh/ehh_scans/chr",chr,".txt"), row.names = FALSE, quote = FALSE)
  
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

#load(file = "data/rehh/test.ehh.scan.Rdata")
test.ihs <- ihh2ihs(test.ehh.scan, freqbin=0.025)

calc_ehhs()

ihsplot(test.ihs$res.ihs, plot.pval=TRUE, ylim.scan=2, main="")

ihs.df <- data.frame(test.ihs$res.ihs)

ihs.df %>% filter(iHS < -3)
(test.in@position == 6771135) %>% which

distribplot(test.ihs$res.ihs[,3])

plot(ihs.df$POSITION, ihs.df$iHS, col = (ihs.df$Pvalue < 0.01)+1)

bifurcation.diagram(test.in, mrk_foc = 680, all_foc = 1, nmrk_l = 2, nmrk_r = 50,
                    main="")


ies2rsb(hh_pop1, hh_pop2, popname1=NA,popname2=NA,method="bilateral")

plot(wht_cmn_rsb[,2], wht_cmn_rsb[,3], col = (wht_cmn_rsb[,4] > 2) +1)

wht_cmn_rsb %>%
  filter(`Pvalue (bilateral)` < -3 | `Pvalue (bilateral)` > 3)

which(wht_cmn_rsb$pos == 12850192)

rsbplot(wht_cmn_rsb)

res.ehh<-calc_ehh(hap.in ,mrk = 787)

bifurcation.diagram(hap.in, mrk_foc = 1461, all_foc = 1, nmrk_l = 5, nmrk_r = 100,
                    main="")
