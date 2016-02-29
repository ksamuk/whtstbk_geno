compute_lc_distance_stats_file <- function (row_num, fst.dat, bathy1, trans1, pop.dat, plot.map = TRUE){
  
  fst_row <- fst.dat[row_num,] %>%
    select(pop1, pop2)
  
  
  pop1 <- fst_row[,1]
  pop2 <- fst_row[,2]
  
	# find coordinates
	
	pop1.coord <- pop.dat %>%
		filter(pop == pop1) %>%
		select(lat, long) %>%
		unique %>%
		unlist %>% 
		as.numeric
	pop1.lat <- pop1.coord[1]
	pop1.lon <- pop1.coord[2]
	
	pop2.coord <- pop.dat %>%
	  filter(pop == pop2) %>%
	  select(lat, long) %>%
	  unique %>%
	  unlist %>% 
	  as.numeric
	pop2.lat <- pop2.coord[1]
	pop2.lon <- pop2.coord[2]
	
	# import NOAA coast data
		
		loc <- data.frame(x = c(pop1.lon, pop2.lon), y = c(pop1.lat, pop2.lat))
		
		# find point to nearest coastline for both pops
		nearest.coastline <- point_to_nearest_coastline(bathy1, loc)
		
		loc$x[1] <- nearest.coastline$loc.x.1.[1] 
		loc$x[2] <- nearest.coastline$loc.x.2.[1]
		loc$y[1] <- nearest.coastline$loc.y.1.[1]
		loc$y[2] <- nearest.coastline$loc.y.2.[1] 
		dist.to.coast1 <- nearest.coastline$dist.to.coast1
		dist.to.coast2 <- nearest.coastline$dist.to.coast2
		
	
	# find the least cost path between coastline points
	least.cost.path <- NA
	least.cost.distance <- NA
	
	try(least.cost.path <- lc_dist_no_bar(trans1, loc), silent = TRUE)
	try(least.cost.distance <- lc_path_to_dist(least.cost.path), silent = TRUE)
	
	# find the great circle distances between the *original* points
	loc.orig <- data.frame( x = c(pop1.lon, pop2.lon), y = c(pop1.lat, pop2.lat)) 
	euc.distance <- earth.dist(loc.orig)[1] %>% round
	
	if (plot.map){
		
		# plot bathy map with lc path + points
		png(paste0("metadata/map_out/",pop1,"_",pop2,".map.png"))
		plot(bathy1, image = TRUE, asp = 1, 
				 land = TRUE, deep= -100000, 
				 shallow=-100, step=100000, 
				 drawlabels = FALSE, 
				 bpal = list(c(min(bathy1,na.rm = TRUE), 0, blues), c(0, max(bathy1, na.rm = TRUE), greys)), 
				 lwd = 0.0)
		
		if (is.na(least.cost.path)){
			dummy <- lapply(loc, lines, col = col2alpha("orange", 0.5), lwd =5, lty = 1)
		}else{
			dummy <- lapply(least.cost.path, lines, col = col2alpha("orange", 0.5), lwd =5, lty = 1)
		}
		points(loc, bg = "orange", cex = 1, pch = 16)
		dev.off()
	}
	
	
	row.out <- data.frame(pop1, pop2, least.cost.distance, euc.distance, dist.to.coast1, dist.to.coast2) %>%
	  distinct
	
	
	
}