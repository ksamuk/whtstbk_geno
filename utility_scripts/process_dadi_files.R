# parse a dadi result output file (of my own invention)

# libraries

library("dplyr")
library("stringr")
library("tidyr")
library("ggplot2")

# Ks estimated by Baocheng Guo in BMC Genomics 2013
mu1 <- 7.1*(10)^-9
mu2 <- 6.8*(10)^-8
theta <- 9723.27623947
Nref <-  theta / (4*mu1*30000*10000*0.8)
Nref <-  theta / (4*mu1*20000*10000*0.7)
div_t <- (T1*2*Nref * 1)/1e6
Nref = theta/(4*mu*L)

# list of dadi results

dadi_files <- list.files("data/dadi/results", full.names = TRUE)

parse_dadi_file <- function(x){
  
  #open file and read lines
  con <- file(x, "r")
  dadi_lines <- readLines(con, warn = FALSE)
  close(con)
  
  pop <- strsplit(x, split = "_") %>% unlist %>% .[3]
  run <- strsplit(x, split = "_") %>% unlist %>% .[4]
  
  #fix weird line break
  dadi_lines[6] <- paste0(dadi_lines[6],dadi_lines[7])
  
  pattern <- "([^0-9. \\-])"
  
  # extract model run info
  grid <- gsub(pattern ,"", dadi_lines[2]) %>% str_trim()
  upper <- gsub(pattern ,"", dadi_lines[3]) %>% str_trim()
  lower <- gsub(pattern ,"", dadi_lines[4]) %>% str_trim()
  initials <- gsub(pattern ,"", dadi_lines[5]) %>% str_trim()
  best_fit <- gsub(pattern ,"", dadi_lines[6]) %>% gsub("^- ", "", .) %>% str_trim()
  log_lik <- gsub(pattern ,"", dadi_lines[8]) %>% str_trim()
  theta <- gsub(pattern ,"", dadi_lines[9]) %>% str_trim()
  
  # format best fit 
  # goes like: 
  # s: intital size pop1 (pop2 = anc - pop1), 
  # nu1: final pop1 size
  # nu2: final pop2 size 
  # T: divergence time 
  # m12: migration from pop1 to pop2
  # m21: migration from pop2 to pop1
  best_fit <- strsplit(best_fit, split = " ") %>% unlist %>% as.character %>% as.numeric %>% na.omit
  
  dadi_df <- data.frame(pop, run, s = best_fit[1], nu1 = best_fit[2], nu2 = best_fit[3],
                        div_t = best_fit[4], m12 = best_fit[5], m21 = best_fit[6], 
                        log_lik, theta, grid, upper, lower, initials)
  dadi_df
}


dadi_df <- lapply(dadi_files, parse_dadi_file)
dadi_df <- bind_rows(dadi_df )

head(dadi_df)

dadi_ln <- dadi_df %>%
  select(-initials, -lower, -upper, -grid, -theta) %>%
  gather(key = coefficient, value = value, -log_lik, -pop, -run)

dadi_ln %>%
  filter(pop == "CL") %>% 
  ggplot(aes (x = as.numeric(run), y = value, color = as.numeric(log_lik)))+
  geom_point()+
  geom_line()+
  facet_grid(coefficient~., scales = "free_y")


dadi_ln %>%
  ggplot(aes (x = pop, y = value))+
  geom_jitter(width = 0.1)+
  #geom_line()+
  facet_grid(coefficient~., scales = "free")


facet_grid(coefficient~., scales = "free")

dadi_ln %>%
  ggplot(aes (x = as.numeric(run), y = value, color = coefficient))+
  geom_point(size = 2)+
  geom_line()+
  facet_grid(coefficient~pop, scales = "free")

