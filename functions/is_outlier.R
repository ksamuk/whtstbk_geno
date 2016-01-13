### identify 95th percentile outliers 

is.outlier <- function(x){

		x95 <- quantile(x, na.rm = TRUE, probs = 0.95)[1]
		return(x >= x95)
}

