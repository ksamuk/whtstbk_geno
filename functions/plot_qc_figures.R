plot_qc_figures <- function(info_df, frac = 0.25){
	
	alpha <- 1 - frac*2
	if (alpha < 0){
		alpha <- 0.1
	}
	
	dp_plot <- info_df %>% 
		sample_frac(frac) %>%
		ggplot(aes(x = pos, y = DP)) +
		geom_point(alpha = alpha, size = 0.5)+
		stat_smooth(method = "loess", se = FALSE, n = 200, span = 0.1, color = "red")+
		scale_y_log10()+
	  xlab("")
	
	dp_hist <- info_df %>% 
		sample_frac(frac) %>%
		ggplot(aes(x = DP)) +
		geom_histogram(bins = 50) +
		scale_x_log10()+
	  ylab("")
	  
	
	mq_plot <- info_df %>% 
		sample_frac(frac) %>%
		ggplot(aes(x = pos, y = MQ)) +
		geom_point(alpha = alpha, size = 0.5)+
		stat_smooth(method = "loess", se = FALSE, n = 200, span = 0.1, color = "red")+
	  xlab("")
	
	mq_hist <- info_df %>% 
		sample_frac(frac) %>%
		ggplot(aes(x = MQ)) +
		geom_histogram(bins = 50)+
	  ylab("")
	
	qd_plot <- info_df %>% 
		sample_frac(frac) %>%
		ggplot(aes(x = pos, y = QD)) +
		geom_point(alpha = alpha, size = 0.5)+
		stat_smooth(method = "loess", se = FALSE, n = 200, span = 0.1, color = "red")+
	  xlab("")
	
	qd_hist <- info_df %>% 
		sample_frac(frac) %>%
		ggplot(aes(x = QD)) +
		geom_histogram(bins = 50)+
	  ylab("")
	
	ind_plot <- info_df %>% 
		mutate(PI = AN/max(AN)) %>% 
		sample_frac(frac) %>%
		ggplot(aes(x = pos, y = PI)) +
		geom_point(alpha = alpha, size = 0.5) +
		stat_smooth(method = "loess", se = FALSE, n = 200, span = 0.1, color = "red")+
	  xlab("")
	
	ind_hist <- info_df %>% 
		mutate(PI = AN/max(AN)) %>% 
		sample_frac(frac) %>%
		ggplot(aes(x = PI)) +
		geom_histogram(bins = 50)+
	  ylab("")
	
	fis_plot <- info_df %>% 
		sample_frac(frac) %>%
		ggplot(aes(x = pos, y = InbreedingCoeff)) +
		geom_point(alpha = alpha, size = 0.5) +
		stat_smooth(method = "loess", se = FALSE, n = 200, span = 0.1, color = "red")+
	  xlab("")
	
	fis_hist <- info_df %>% 
		sample_frac(frac) %>%
		ggplot(aes(x = InbreedingCoeff)) +
		geom_histogram(bins = 50)+
	  ylab("")
	
	plot_grid(dp_plot, dp_hist, 
						mq_plot, mq_hist, 
						qd_plot, qd_hist, 
						ind_plot, ind_hist,
						fis_plot, fis_hist,
						align = "hv", ncol = 2)
	
}