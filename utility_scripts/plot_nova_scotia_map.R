library(raster)
library(rgdal)
library(ggplot2)
library(reshape2)
library(plyr)

#  Assuming you have a path 'Maps' that you store your spatial files in.  This
#  is all downloaded from <a href="http://www.naturalearthdata.com/downloads/">http://www.naturalearthdata.com/downloads/</a> using the
#  1:50m "Medium" scale data.

nat.earth <- stack('data/maps/NE1_HR_LC_SR_W_DR/NE1_HR_LC_SR_W_DR.tif')

ne_coast <- readOGR('data/maps/NE1_HR_LC_SR_W_DR/ne_10m_coastline.shp','ne_10m_coastline')

ne_rivers <- readOGR('./Maps/NaturalEarth/ne_50m_rivers_lake_centerlines.shp',
                     'ne_50m_rivers_lake_centerlines')

ne_coast <- readOGR('./Maps/NaturalEarth/ne_50m_coastline.shp',
                    'ne_50m_coastline')

#  I have a domain I'm interested in, but there's no reason you can't define something else:
quick.subset <- function(x, longlat){
  
  # longlat should be a vector of four values: c(xmin, xmax, ymin, ymax)
  x@data$id <- rownames(x@data)
  
  x.f = fortify(x, region="id")
  x.join = join(x.f, x@data, by="id")
  
  x.subset <- subset(x.join, x.join$long > longlat[1] & x.join$long < longlat[2] &
                       x.join$lat > longlat[3] & x.join$lat < longlat[4])
  
  x.subset
}


#43°21'53.9"N 66°44'27.4"W
#47°12'31.0"N 60°03'24.8"W
#domain <- c(-98.6, -66.1, 36.5, 49.7)
domain <- c(-66.44,-60.44, 42.21, 47.12)
#domain <- c(42.21, -66.44, 47.12, -60.44)
#domain <- c(42.21, -66.44, 47.12, -60.44)
coast.subset <- quick.subset(ne_coast, domain)

nat.crop <- crop(nat.earth, y=extent(domain))

rast.table <- data.frame(xyFromCell(nat.crop, 1:ncell(nat.crop)),
                         getValues(nat.crop/255))

rast.table$rgb <- with(rast.table, rgb(NE1_HR_LC_SR_W_DR.1,
                                       NE1_HR_LC_SR_W_DR.2,
                                       NE1_HR_LC_SR_W_DR.3,
                                       1))
# et voila!

ggplot(data = rast.table, aes(x = x, y = y)) +
  geom_tile(fill = rast.table$rgb) +
  #geom_polygon(data=lakes.subset, aes(x = long, y = lat, group = group), fill = '#ADD8E6') +
  scale_alpha_discrete(range=c(1,0)) +
  #geom_path(data=lakes.subset, aes(x = long, y = lat, group = group), color = 'blue') +
  #geom_path(data=river.subset, aes(x = long, y = lat, group = group), color = 'blue') +
  geom_path(data=coast.subset, aes(x = long, y = lat, group = group), color = 'lightblue') +
  scale_x_continuous(expand=c(0,0)) +
  scale_y_continuous(expand=c(0,0)) +
  xlab('') + ylab('')