library(rgdal)
library(RColorBrewer)
library("viridis")
dpath<-"C:/Users/Kieran/Desktop/ns maps/e055ns20/grid/dem020hy/w001001.adf"
x <- new("GDALReadOnlyDataset", dpath)
getDriver(x)
getDriverLongName(getDriver(x))


grad <- seq(0, 550, by = 10)
pal <- viridis(length(grad))


xx<-asSGDF_GROD(x, output.dim=c(4000, 4000))
spplot(xx, "band1", at=grad, col.regions=pal)
