# classify longitudes into regions
# based on stickleback populations from Samuk et al. 2015, so possibly regions are:
# NA, EU, JP

lon_to_region <- function(lon){
	
	if(lon < -20){
		region <- "na"
	} else if(lon > 100){
		region <- "jp"
	} else if(lon < 50 & lon > -20){
		region <- "eu"
	}
	
	return(region)
	
}