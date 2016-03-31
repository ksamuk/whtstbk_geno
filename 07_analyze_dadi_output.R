# parse a dadi result output file (of my own invention)

################################################################################
# initials
################################################################################

library("dplyr")
library("stringr")
library("tidyr")
library("ggplot2")
library("ggthemes")

################################################################################
# initials
################################################################################

# list of dadi results

dadi_files <- list.files("data/dadi/results_final", full.names = TRUE)


################################################################################
# parse dadi function 
################################################################################

parse_dadi_file <- function(x, mode = "im"){
  
  if(mode != "split"){
    #open file and read lines
    con <- file(x, "r")
    dadi_lines <- readLines(con, warn = FALSE)
    close(con)
    
    pop <- strsplit(x, split = "_") %>% unlist %>% .[4]
    run <- strsplit(x, split = "_") %>% unlist %>% .[5]
    
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
    
    dadi_df <- data.frame(pop, run, model = "im", s = best_fit[1], nu1 = best_fit[2], nu2 = best_fit[3],
                          div_t = best_fit[4], m12 = best_fit[5], m21 = best_fit[6], 
                          log_lik, theta, grid, upper, lower, initials)
    dadi_df
  }else{
    #open file and read lines
    con <- file(x, "r")
    dadi_lines <- readLines(con, warn = FALSE)
    close(con)
    
    pop <- strsplit(x, split = "_") %>% unlist %>% .[4]
    run <- strsplit(x, split = "_") %>% unlist %>% .[6] %>% gsub("mig", "", .)
    
    #fix weird line break
    #dadi_lines[6] <- paste0(dadi_lines[6],dadi_lines[7])
    
    pattern <- "([^0-9. \\-])"
    
    # extract model run info
    grid <- gsub(pattern ,"", dadi_lines[2]) %>% str_trim()
    upper <- gsub(pattern ,"", dadi_lines[3]) %>% str_trim()
    lower <- gsub(pattern ,"", dadi_lines[4]) %>% str_trim()
    initials <- gsub(pattern ,"", dadi_lines[5]) %>% str_trim()
    best_fit <- gsub(pattern ,"", dadi_lines[6]) %>% gsub("^- ", "", .) %>% str_trim()
    log_lik <- gsub(pattern ,"", dadi_lines[7]) %>% str_trim()
    theta <- gsub(pattern ,"", dadi_lines[8]) %>% str_trim()
    
    # format best fit 
    # goes like: 
    # s: intital size pop1 (pop2 = anc - pop1), 
    # nu1: final pop1 size
    # nu2: final pop2 size 
    # T: divergence time 
    # m12: migration from pop1 to pop2
    # m21: migration from pop2 to pop1
    best_fit <- strsplit(best_fit, split = " ") %>% unlist %>% as.character %>% as.numeric %>% na.omit
    
    dadi_df <- data.frame(pop, run, model = "split_mig", s1 = best_fit[1], s2 = best_fit[2], div_t = best_fit[3],
                          m = best_fit[4], log_lik, theta, grid, upper, lower, initials)
    dadi_df
  }
  
  
}

################################################################################
# process dadi files into data frames
################################################################################

# parse all files and bind to data frame

dadi_files_im <- grep("split", dadi_files, invert = TRUE, value = TRUE)
dadi_files_split <- grep("split", dadi_files, value = TRUE)

dadi_df_im <- lapply(dadi_files_im, parse_dadi_file)
dadi_df_im <- bind_rows(dadi_df_im)

dadi_df_split <- lapply(dadi_files_split, parse_dadi_file, mode = "split")
dadi_df_split <- bind_rows(dadi_df_split)

dadi_im_ln <- dadi_df_im  %>%
  select(-initials, -lower, -upper, -grid, -theta) %>%
  gather(key = coefficient, value = value, -log_lik, -pop, -run, -model)

dadi_split_ln <- dadi_df_split  %>%
  select(-initials, -lower, -upper, -grid, -theta) %>%
  gather(key = coefficient, value = value, -log_lik, -pop, -run, -model)

################################################################################
# process dadi files into data frames
################################################################################

# what are the best fit im models for each pair of populations?

dadi_df_im %>% 
  mutate(log_lik = as.numeric(log_lik)) %>%
  filter(pop %in% c("CL", "SR")) %>%
  #filter(pop %in% c("SR")) %>%
  filter(!is.na(m21)) %>%
  group_by(pop) %>%
  filter(log_lik == max(log_lik))

#pop   run  model         s       nu1       nu2     div_t      m12      m21   log_lik         theta     grid
#(chr) (chr) (fctr)     (dbl)     (dbl)     (dbl)     (dbl)    (dbl)    (dbl)     (dbl)         (chr)    (chr)
#1    CL    21     im 0.2757301 0.1468878 0.2133241 0.5780141 39.61412 21.10621 -1980.307 10521.1294381 50 60 70
#2    SR   18a     im 0.1744412 0.2025371 0.3132589 1.1991450 39.60331 18.91881 -2641.557 10048.3453029 50 60 70

# what are the best fit 'split-mig' models for each pair of populations?

dadi_df_split %>% 
  mutate(log_lik = as.numeric(log_lik)) %>%
  filter(pop %in% c("CL", "SR")) %>%
  group_by(pop) %>%
  filter(log_lik == max(log_lik))

# Source: local data frame [2 x 13]
# Groups: pop [2]
# 
# pop   run     model        s1        s2    div_t        m   log_lik         theta     grid         upper
# (chr) (chr)    (fctr)     (dbl)     (dbl)    (dbl)    (dbl)     (dbl)         (chr)    (chr)        (fctr)
# 1    CL    16 split_mig 0.4463100 0.3415070 16.06134 6.484532 -5940.279 8020.08576759 50 60 70 0.5 0.5 50 20
# 2    SR    18 split_mig 0.3909826 0.4558576 14.45678 6.544098 -6114.416  7471.2396129 50 60 70 0.5 0.5 50 20

# CL SNM VALUES:

#Maximum log composite likelihood: -3765.07988701
#Optimal value of theta: 6536.83971583

# SR SNM VALUES:

#Maximum log composite likelihood: -3571.07453008
#Optimal value of theta: 6577.23873737

# note on units:
# Times are given in units of 2Nref generations
# Migration rates are given in units of Mij = 2Nrefmij

################################################################################
# perform the likelihood ratio tests
################################################################################

# CL
log_lik_im <-  -1980.307
log_lik_split_mig <-  -5940.279 
log_lik_snm <-  -3765.07988701

# im vs. snm
# im = 7 params, snm = 1 param (theta)
df <- 6 
chi_sq <- (log_lik_im - log_lik_snm)*2
pchisq(chi_sq, df = df, lower.tail = FALSE)  %>% as.numeric()

# im vs. splitmig
# im = 7 params, splut = 5 param 
df <- 2 
chi_sq <- (log_lik_im - log_lik_split_mig)*2
pchisq(chi_sq, df = df, lower.tail = FALSE)  %>% as.numeric()

# SR
log_lik_im <- -3325.019
log_lik_split_mig <- -6114.416
log_lik_snm <- -3571.07453008

# im vs. snm
# im = 7 params, snm = 1 param (theta)
df <- 6 
chi_sq <- (log_lik_im - log_lik_snm)*2
pchisq(chi_sq, df = df, lower.tail = FALSE)  %>% as.numeric()

# im vs. splitmig
# im = 7 params, splut = 5 param 
df <- 2 
chi_sq <- (log_lik_im - log_lik_split_mig)*2
pchisq(chi_sq, df = df, lower.tail = FALSE)  %>% as.numeric()

################################################################################
# transformed parameters 
################################################################################

# pop   run  model         s       nu1       nu2      div_t      m12      m21   log_lik         theta
# (chr) (chr) (fctr)     (dbl)     (dbl)     (dbl)      (dbl)    (dbl)    (dbl)     (dbl)         (chr)
# 1    CL    21     im 0.2757301 0.1468878 0.2133241  0.5780141 39.61412 21.10621 -1980.307 10521.1294381
# 2    SR    28     im 0.4416730 0.5793114 0.5266125 11.4618391 19.46700 22.06862 -3325.019 5876.76388284


# Ks estimated by Baocheng Guo in BMC Genomics 2013
mu1 <- 7.1*(10)^-9
# mu estimated by Marius Roesti in Nat Comm. 2016
mu2 <- 6.8*(10)^-8

# total number of sequenced sites (not snps) *after* max-missing filter
effective_length <- 530125

# CL estimates

theta <- 10521.1294381
t1 <- 0.5780141
m12 <- 39.61412
m21 <- 21.10621
downsampling <- 0.81 # proportion of sites filtered/lost during projection/etc.

n_ref <-  theta / (4*mu1*effective_length*downsampling)

# "Times are given in 2Nref generations"
# If Nref=20000(individuals):
# Real time =2*Nref*T    , suppose T is inferred by dadi and T equals 1, 
# and the generation time is 1 year. Then Real time equals 40000 years.

div_t <- (t1*2*n_ref * 1)

# "Migration rates are given in 2Nref*migration rate"
# if nref = 20000
# "mij does equal 2.5e-4, but that is not the number of individuals per year. 
# Instead it is the fraction of individuals in each generation in population i 
# who are new migrants from population j.
# So if the size of population i is 0.5 Nref, then 2.5 individuals each generation are migrating from population j to population i."
nu1 <- 0.1468878
nu2 <- 0.2133241

pop_size1 <- nu1*n_ref
pop_size2 <- nu2*n_ref

mig_1_2 <- (m12/(2*n_ref))*pop_size1
mig_2_1 <- (m21/(2*n_ref))*pop_size2 

################################################################################
# plots of model uncertainty
################################################################################

dadi_im_ln %>%
  filter(pop %in% c("CL", "SR")) %>%
  ggplot(aes (x = pop, y = value))+
  geom_jitter(width = 0.1, size = 1)+
  stat_summary(fun.y=median, fun.ymin=median, fun.ymax=median, 
               geom="crossbar", width=0.3, color = "red", fatten = 2)+
  facet_wrap(~coefficient, scales = "free")+
  ylim(0, NA)+
  theme_base() +
  theme(legend.background = element_blank(),
        text = element_text(size = 14),
        axis.title.x = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        legend.justification = c(1, 1), 
        legend.position = "none",
        #axis.title.y = element_text(margin=margin(0,6,0,0), face = "bold"),
        #axis.title.x = element_text(margin=margin(12,0,0,0), face = "bold"),
        plot.background = element_rect(color = NA))+
  scale_fill_fivethirtyeight()


dadi_split_ln %>%
  filter(pop %in% c("CL", "SR")) %>%
  ggplot(aes (x = pop, y = value))+
  geom_jitter(width = 0.1)+
  #geom_line()+
  facet_grid(coefficient~., scales = "free")+
  ylim(0,NA)

dadi_im_ln %>%
  filter(pop != "BS") %>%
  ggplot(aes (x = pop, y = value))+
  geom_jitter(width = 0.1)+
  #geom_line()+
  facet_grid(coefficient~., scales = "free")+
  ylim(0,NA)


facet_grid(coefficient~., scales = "free")

dadi_ln %>%
  ggplot(aes (x = as.numeric(run), y = value, color = coefficient))+
  geom_point(size = 2)+
  geom_line()+
  facet_grid(coefficient~pop, scales = "free")

################################################################################
# converting model results to biologically relevant units
################################################################################

# Ks estimated by Baocheng Guo in BMC Genomics 2013
mu1 <- 7.1*(10)^-9
# mu estimated by Marius Roesti in Nat Comm. 2016
mu2 <- 6.8*(10)^-8

# total number of sequenced sites (not snps) *after* max-missing filter
effective_length <- 530125

# effective length = 530125 (determined by emmitting all including invariant + counting called genotypes)

theta <- 6577
Nref <-  theta / (4*mu1*530125*0.51)

# "2Nref generations"
div_t <- (T1*2*Nref * 1)/1e6
Nref = theta/(4*mu*L)

