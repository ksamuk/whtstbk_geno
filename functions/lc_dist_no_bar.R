# adaptation of lc.dist from 'marmap'
# no progress bar, no dist option

lc_dist_no_bar <- function (trans, loc) {

		nb.loc <- nrow(loc)
		path <- list()
		comb <- combn(1:nb.loc, 2)
		for (i in 1:ncol(comb)) {
			origin <- sp::SpatialPoints(loc[comb[1, i], ])
			goal <- sp::SpatialPoints(loc[comb[2, i], ])
			temp <- gdistance::shortestPath(trans, origin, goal, 
																			output = "SpatialLines")
			path[[i]] <- temp@lines[[1]]@Lines[[1]]@coords
		}
		return(path)
	}
