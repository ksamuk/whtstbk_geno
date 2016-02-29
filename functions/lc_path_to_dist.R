# calculate total distance from an lc.dist object

lc_path_to_dist <- function(least.cost.path){
	path <- least.cost.path[[1]]
	dist.vec <- rep(NA, length(path[,1]))
	
	for (i in 2:length(path[,1])){
		
		long1 <- path[,1][i]
		lat1 <- path[,2][i]
		long2 <- path[,1][i-1]
		lat2 <- path[,2][i-1]
		
		dist.vec[i] <- deg.dist(long1, lat1, long2, lat2)
		
	}
	
	return(sum(dist.vec, na.rm = TRUE))
}