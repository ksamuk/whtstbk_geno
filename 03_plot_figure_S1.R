rm(list=ls())

################################################################################
# Libraries and initalizing variables
################################################################################

library("mapdata")
library("maps")
library("maptools")
library("rworldmap")
library("dplyr")
library("mapplots")
library("plotrix")
library("ggthemes")
library("ggplot2")
library("gridBase")
library("grid")
library("tidyr")

list.files("shared_functions", full.names = TRUE) %>% sapply(source) %>% invisible

################################################################################
# Plot the map
################################################################################

pop.dat <- read.csv("metadata/whtstbk_site_coordinates.csv", header = TRUE, stringsAsFactors = FALSE)
margins_maps <- c(2,2,2,2)

pdf(file = "figures/Figure1.pdf", width = 8.5, height = 8.5)
par(mfrow = c(1,2), mar = c(1,1,1,1))

# whole provence 
map("worldHires","Canada", xlim=c(-66.5,-59.75), ylim=c(43.25,47.25), col="gray90", fill=TRUE,
    resolution = 0, mar = margins_maps)
map.axes()

lapply(1:length(pop.coords[[1]]), plot_pies, pop.coords = pop.coords, radius = 0.12, cex = 1.8) %>%  invisible()
lapply(1:length(pop.coords[[1]]), plot_labels,  pop.dat =  pop.dat) %>%  invisible()

rect(xleft = -62.125, ybottom = 45.25, xright = -60.5, ytop = 46.2)

map.scale(x = -62, y= 44, ratio = FALSE, cex = 0.7 )


# straight of canso region
map("worldHires","Canada", xlim=c(-62.125,-60.5), ylim=c(45.25,46.2), col="gray90", fill=TRUE,
    resolution = 0, mar = margins_maps)
map.axes()

lapply(1:length(pop.coords[[1]]), plot_pies, pop.coords = pop.coords, radius = 0.035, cex = 1.9, zoom = TRUE) %>%  invisible()
lapply(1:length(pop.coords[[1]]), plot_labels,  pop.dat =  pop.dat, zoom = TRUE) %>%  invisible()

map.scale(x = -60.85, y= 45.4, ratio = FALSE, cex = 0.7 )

dev.off()