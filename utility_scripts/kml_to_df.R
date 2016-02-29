# convert a kml file to geographic coordinates

#install.packages("maptools")
library("maptools")
library("rgdal")
library("sp")
library("dplyr")

######### READ KML FUNCTION
# from: https://gist.github.com/holstius/6631918
read.kml <- function(file, layers) {
  require(sp)
  require(rgdal)
  read.layer <- function (layer_name) {
    spobj <- rgdal::readOGR(dsn=file, layer=layer_name)
    coords <- coordinates(spobj)
    colnames(coords) <- c('x', 'y', 'z')[1:ncol(coords)]
    df <- data.frame(coords, spobj@data)
    transform(df, layer=layer_name)
  }
  Reduce(rbind, lapply(layers, read.layer))
}
#########


kml_file <- "metadata/whtstbk_geno_collection_sites.kml"

kml_df <- read.kml(file = kml_file, layers = "Final Sites") %>%
  select(-Description, -layer, -z)

write.csv(kml_df, "metadata/whtstbk_site_coordinates.csv")


