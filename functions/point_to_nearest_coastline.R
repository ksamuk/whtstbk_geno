
# move a point to the nearest coastline

point_to_nearest_coastline <- function (bat, loc) {
  nearest.coastline <- NA
  dist.to.coast1 <- NA
  dist.to.coast2 <- NA
  
  try(nearest.coastline <- dist2isobath(bat, loc, isobath = 10))
  
  if (!is.na(nearest.coastline[, 1][1])) {
    loc$x <- nearest.coastline[, 4]
    loc$y <- nearest.coastline[, 5]
    dist.to.coast1 <- nearest.coastline[, 1][1] / 1000 # meters
    dist.to.coast2 <- nearest.coastline[, 1][2] / 1000 # meters
    
    
  }
  
  return(
    data.frame(
      nearest.coastline,
      dist.to.coast1,
      dist.to.coast2,
      loc$x[1],
      loc$y[1],
      loc$x[2],
      loc$y[2]
    )
  )
}
