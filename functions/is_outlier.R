### identify 95th percentile outliers 

is.outlier <- function(x, cutoff = 0.95){

		x95 <- quantile(x, na.rm = TRUE, probs = cutoff)[1]
		return(x >= x95)
}

