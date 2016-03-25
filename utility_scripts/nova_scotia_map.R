

#min     max
#x  210068  767668library(rgdal)
library(raster)
library(RColorBrewer)
library("viridis")
library("dplyr")
dpath<-"C:/Users/Kieran/Desktop/ns maps/e055ns20/grid/dem020hy/w001001.adf"
x <- new("GDALReadOnlyDataset", dpath)
getDriver(x)
getDriverLongName(getDriver(x))


grad <- seq(1, 550, by = 2)
pal <- viridis(length(grad))
pal <- colorRampPalette(c("black","white"),length(grad))
#y 4797718 5239478

xx<-asSGDF_GROD(x, output.dim=c(5000, 5000))
tmp <- spplot(xx, "band1", at=grad, col.regions=pal, xlim = c(520000, 750000), ylim =c(4999478,5139478), axes=FALSE)
png("myplot.png", width = 5000, height = 5000)
tmp
dev.off()

bright_tiff <- "data/bright_atlas.geotiff"
sat <- raster(bright_tiff)
plot(sat,xmin = -66, xmax = -60, ymin = 140, ymax = 144)
