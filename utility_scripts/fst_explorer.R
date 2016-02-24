

pop1 <- c("AAAT", "AAAT", "AAAT", "AAAT", "AAAT", "AAAA", "AAAA", "AAAA")
pop2 <- c("AAAA", "AAAA", "AAAA", "AAAA", "AAAA", "AAAA", "AAAA", "AAAT")

calc_fst(pop1, pop2)
calc_dxy(pop1, pop2)

calc_dxy <- function(pop1, pop2){
  sapply(pop1, function(x) adist(x, pop2) %>% mean) %>% mean
}


calc_fst <- function(pop1, pop2){
  
  fst <- list()
  
  for (i in 1:nchar(pop1[1])){
    
    # extract pop1 alleles
    pop1_alleles <- lapply(pop1, function(x)substr(x, start = i, stop = i)) %>% unlist
    
    # extract pop2 alleles
    pop2_alleles <- lapply(pop2, function(x)substr(x, start = i, stop = i)) %>% unlist
    
    # all alleles
    all_alleles <- c(pop1_alleles, pop2_alleles)
    
    # ht, aka 2pq
    
    # if site is monomorphic, fst = 0
    if (length(unique(all_alleles)) == 1){
      fst[[i]] <- data.frame(fst = Inf, ht = 0, hs = 0)
    } else{
      
      # count the alleles in total population
      allele_counts <- table(all_alleles) %>% data.frame
      
      # calc frequencies (total)
      p <- allele_counts$Freq[1] / sum(allele_counts$Freq)
      q <- allele_counts$Freq[2] / sum(allele_counts$Freq)
      
      # overall expected heterozygosity
      ht <- 2*p*q
      
      # average intrapopulation heterozygostity
      pop1_counts <- table(pop1_alleles) %>% data.frame
      
      if (nrow(pop1_counts) == 1){
        hs_pop1 <- 0
      } else{
        p_pop1 <-  pop1_counts$Freq[1] / sum(pop1_counts$Freq)
        q_pop1 <-  pop1_counts$Freq[2] / sum(pop1_counts$Freq)
        
        hs_pop1 <- 2*p_pop1*q_pop1
      }
      
      pop2_counts <- table(pop2_alleles) %>% data.frame
      
      if (nrow(pop2_counts) == 1){
        hs_pop2 <- 0
      } else{
        p_pop2 <-  pop2_counts$Freq[1] / sum(pop2_counts$Freq)
        q_pop2 <-  pop2_counts$Freq[2] / sum(pop2_counts$Freq)
        
        hs_pop2 <- 2*p_pop2*q_pop2
      }
      
      hs <- mean(c(hs_pop1, hs_pop2))
      
      # fst!
      fst[[i]] <- data.frame(fst=(ht-hs)/ht, ht = ht, hs = hs)
      
    }
    
  }
  return(fst)
  
}