rm(list=ls())

################################################################################
# Libraries and initalizing variables
################################################################################

library("mapdata")
library("maps")
library("maptools")
library("rworldmap")
library("dplyr")
library("plotrix")

list.files("shared_functions", full.names = TRUE) %>% sapply(source) %>% invisible

################################################################################
# Plot figure 1
################################################################################

# input files
# locations of nova scotia populations
pop.dat <- read.csv("metadata/whtstbk_site_coordinates.csv", header = TRUE, stringsAsFactors = FALSE)


# population points
# top and bottom are for the half-circle points
pop.coords <- list(x = pop.dat$long, y = pop.dat$lat)
pop.coords.bottom <- list(x = pop.dat$long, y = pop.dat$lat-0.025)

################################################################################
# Pie plotting function
################################################################################

plot_pies <- function(index, pop.coords, radius = 0.1, cex = 2.3, zoom = FALSE){
  
  if(pop.dat$ecotypes[index] == "both"){
    halves <- c(0.5,0.5)
    col <- c("white", "blue")
    
    floating.pie(pop.coords$x[index],pop.coords$y[index],halves,
                 edges=500,radius=radius,startpos=1.57,
                 col = col)
    
  } else if (pop.dat$ecotypes[index] == "common"){
    col <- c("blue")
    points(pop.coords$x[index],pop.coords$y[index], pch=21, col = "black", bg = col, cex = cex)

  }else{
    halves <- c(1)
    col <- c("white")
  }
  
}


plot_labels <- function(index, pop.dat, offset_y_default = 0.15, offset_x_default = 0, zoom = FALSE){
  
  if(!pop.dat$plot[index] & zoom == FALSE){
    return("null")
  }
  
  offset_y <- pop.dat$offset_y[index]
  offset_x <- pop.dat$offset_x[index]
  
  offset_y <- ifelse(is.na(offset_y), offset_y_default, offset_y)
  offset_x <- ifelse(is.na(offset_x), offset_x_default, offset_x)
  
  if(zoom == TRUE){
    
    offset_y <- pop.dat$offset_y_zoom[index]
    offset_x <- pop.dat$offset_x_zoom[index]
    
    offset_y <- ifelse(is.na(offset_y), offset_y_default, offset_y)
    offset_x <- ifelse(is.na(offset_x), offset_x_default, offset_x)
    
  }
  
  text(pop.dat$long[index]+offset_x, pop.dat$lat[index]+offset_y, labels = pop.dat$pop[index], font = 2, cex = 0.7)
  
}
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

################################################################################
# Plot the map (cluster colors)
################################################################################

pop.dat <- read.csv("metadata/whtstbk_site_coordinates.csv", header = TRUE, stringsAsFactors = FALSE)
meta_df <- read.csv("metadata/mega_meta.csv", stringsAsFactors = FALSE)

#harmonize names
meta_df$pop[meta_df$pop=="WR"] <- "RR"

# calculate cluster proportions

prop_df <- meta_df %>%
  filter(!is.na(cluster)) %>%
  filter(cluster != "dk") %>%
  group_by(pop, cluster) %>%
  summarise (n = n()) %>%
  mutate(freq = n / sum(n)) %>%
  select(pop, cluster, freq) %>%
  spread(key = cluster, value = freq)

prop_df[is.na(prop_df)] <- 0 

pop.dat <- left_join(pop.dat, prop_df)

################################################################################
# Pie plotting function (clusters)
################################################################################

plot_pies_cluster <- function(index, pop.coords, radius = 0.1, cex = 2.3){
  
  props <- pop.dat[,c("cbr", "cmn", "wht")][index,] %>% as.numeric()

  col <- c("#77AB43","#008FD5","#FFFFFF")
  print(props)
  print(col)
  
  #floating.pie(pop.coords$x[index],pop.coords$y[index], x = props, 
  #             edges=500,radius=radius)
  
  add.pie(z=c(props), x=pop.coords$x[index], y=pop.coords$y[index], radius=radius, 
          col=col, labels="")
  
}

################################################################################
# Plot the map
################################################################################

margins_maps <- c(2,2,2,2)

pdf(file = "figures/FigureS1.pdf", width = 8.5, height = 8.5)
par(mfrow = c(1,2), mar = c(1,1,1,1))

# whole provence 
map("worldHires","Canada", xlim=c(-66.5,-59.75), ylim=c(43.25,47.25), col="gray90", fill=TRUE,
    resolution = 0, mar = margins_maps)
map.axes()

lapply(1:length(pop.coords[[1]]), plot_pies_cluster, pop.coords = pop.coords, radius = 0.1) %>%  invisible()
lapply(1:length(pop.coords[[1]]), plot_labels,  pop.dat =  pop.dat) %>%  invisible()

rect(xleft = -62.125, ybottom = 45.25, xright = -60.5, ytop = 46.2)

map.scale(x = -62, y= 44, ratio = FALSE, cex = 0.7 )


# straight of canso region
map("worldHires","Canada", xlim=c(-62.125,-60.5), ylim=c(45.25,46.2), col="gray90", fill=TRUE,
    resolution = 0, mar = margins_maps)
map.axes()

lapply(1:length(pop.coords[[1]]), plot_pies_cluster, pop.coords = pop.coords, radius = 0.025) %>%  invisible()
lapply(1:length(pop.coords[[1]]), plot_labels,  pop.dat =  pop.dat, zoom = TRUE) %>%  invisible()

map.scale(x = -60.85, y= 45.4, ratio = FALSE, cex = 0.7 )

dev.off()


