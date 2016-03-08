# dat baypass tho?

# libraries
library("dplyr")
library("tidyr")
library("ggplot2")

slug <- "wht_cmn"

list.files("functions", full.names = TRUE) %>% sapply(.,source, verbose = FALSE, echo = FALSE) %>% invisible

# the site labels in their original order
baypass_sites <- read.table("metadata/baypass_sites.txt", header = TRUE, stringsAsFactors = FALSE)

# the baypass output files
xtx_file <- list.files("data/baypass/wht_cmn", pattern = paste0(slug, "_summary_pi_xtx"),full.names = TRUE)
xtx_df <- read.table(xtx_file, header = TRUE, stringsAsFactors = FALSE)
head(xtx_df)

names(xtx_df) <- tolower(names(xtx_df)) 

xtx_df <- cbind(baypass_sites, xtx_df[,-1])

# plot?

xtx_df %>%
  ggplot(aes(x = pos, y = m_xtx, color = as.factor(chr))) +
  geom_point(size = 1 )+
  facet_grid(.~chr, scales = "free", space = "free_x", switch = "x") +
  theme_classic()+
  theme(axis.text.x = element_blank(), 
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        panel.margin = unit(0.1, "cm"),
        strip.background = element_blank(),
        legend.position = "none")


#get estimates (post. mean) of both the a_pi and b_pi parameters of
#the Pi Beta distribution
beta_file <- list.files("data/baypass/wht_cmn", pattern = paste0(slug, "_summary_beta_params"),full.names = TRUE)
pi_beta_coef <- read.table(beta_file, h = TRUE)$Mean

#upload the original data to obtain total allele count
geno_file <- list.files("data/baypass/wht_cmn", pattern = paste0(slug, ".geno"),full.names = TRUE)
geno_data <- geno2YN(geno_file)

# the covariance matrix
omega_file <- list.files("data/baypass/wht_cmn", pattern = paste0(slug, "_mat_omega"),full.names = TRUE)
omega <- as.matrix(read.table(omega_file))


#### RAN POD FILES IN BAYPASS

#######################################################
#Sanity Check: Compare POD and original data estimates
#######################################################

# get estimate of omega from the POD analysis
pod_omega_file <- list.files("data/baypass/wht_cmn", pattern = paste0(slug, "_pods_mat_omega"), full.names = TRUE)
pod_omega <- as.matrix(read.table(pod_omega_file))
plot(pod_omega, omega)
abline(a = 0, b = 1)
fmd.dist(pod_omega, omega)

# get estimates (post. mean) of both the a_pi and b_pi parameters of
# the Pi Beta distribution from the POD analysis
pod_beta_file <- list.files("data/baypass/wht_cmn", pattern = paste0(slug, "_pods_summary_beta_params"),full.names = TRUE)
pod_pi_beta_coef <- read.table(pod_beta_file, h = TRUE)$Mean
plot(pod_pi_beta_coef, pi_beta_coef) ; abline(a = 0, b = 1)

#######################################################
#XtX calibration
#######################################################


pod_xtx_file <- list.files("data/baypass/wht_cmn", pattern = paste0(slug, "_pods_summary_pi_xtx"),full.names = TRUE)
pod_xtx_df <- read.table(xtx_file, header = TRUE, stringsAsFactors = FALSE)
names(pod_xtx_df) <- tolower(names(pod_xtx_df)) 
pod_xtx_df <- cbind(baypass_sites, pod_xtx_df[,-1])

#compute the 1% threshold
pod_thresh <- quantile(pod_xtx_df$m_xtx, probs=0.99)

xtx_df$xtx_outlier <- xtx_df$m_xtx >= pod_thresh

#add the thresh to the actual XtX plot
xtx_df %>%
  ggplot(aes(x = pos, y = m_xtx, color = as.factor(chr))) +
  geom_point(size = 1 )+
  geom_hline(yintercept = pod_thresh) +
  facet_grid(.~chr, scales = "free", space = "free_x", switch = "x") +
  theme_classic()+
  theme(axis.text.x = element_blank(), 
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        panel.margin = unit(0.1, "cm"),
        strip.background = element_blank(),
        legend.position = "none")


#######################################################
#AUX model results (LD control)
#######################################################

aux_beti_file <- list.files("data/baypass/wht_cmn", pattern = paste0(slug, "_aux_summary_betai"),full.names = TRUE)
aux_snp_res <- read.table(aux_beti_file ,h = TRUE)

aux_df <- cbind(baypass_sites, aux_snp_res[,-c(1:2)])

hist(aux_df$BF.dB.)

aux_df %>%
  #filter(BF.dB. > 10) %>%
  ggplot(aes(x = pos, y = BF.dB., color = as.factor(chr))) +
  geom_point(size = 2 )+
  #geom_hline(yintercept = pod_thresh) +
  facet_grid(.~chr, scales = "free_y", space = "free_x", switch = "x") +
  theme_classic()+
  theme(axis.text.x = element_blank(), 
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        panel.margin = unit(0.1, "cm"),
        strip.background = element_blank(),
        legend.position = "none")


graphics.off()
layout(matrix(1:2,2,1))
plot(aux_snp_res$M_Delta,xlab="SNP",ylab=expression(delta[i]),main="AUX model")
plot(covauxisb1.snp.res$M_Delta,xlab="SNP",ylab=expression(delta[i]),
     main="AUX model with isb=1")
