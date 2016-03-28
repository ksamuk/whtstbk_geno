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
# Plot the map (pca + cluster colors)
################################################################################

pop.dat <- read.csv("metadata/whtstbk_site_coordinates.csv", header = TRUE, stringsAsFactors = FALSE)
meta_df <- read.csv("metadata/mega_meta.csv", stringsAsFactors = FALSE, header = TRUE)
pca_df <- read.table("metadata/pca_df.txt", stringsAsFactors = FALSE, header = TRUE)

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
  #print(props)
  #print(col)
  
  #floating.pie(pop.coords$x[index],pop.coords$y[index], x = props, 
  #             edges=500,radius=radius)
  
  add.pie(z=c(props), x=pop.coords$x[index], y=pop.coords$y[index], radius=radius, 
          col=col, labels="")
  
}

################################################################################
# Plot the map
################################################################################

cluster_short <- c("cmn", "cbr", "wht")
cluster_names <- c("Mainland Common", "Bras d'Or Common", "White")

var_prop <- read.table("metadata/pca_var_df.txt", h = T)[,1]

# pca of clusters
cluster_names <- c("Mainland Common", "Bras d'Or Common", "White")
pca_clust <- pca_df %>%
  mutate(Year = as.factor(year) %>% reorder (., . == "2012")) %>%
  mutate(year_pch = c(21,24)[as.numeric(Year)]) %>%
  mutate(Sex = sex) %>%
  mutate(Sex = ifelse(Sex == "M", "Male", "Female")) %>%
  mutate(Group = cluster_names[match(cluster, cluster_short)]) %>%
  #filter(!(id %in% c("SR20", "LN29", "SR15", "2012_SF16"))) %>%
  filter(!is.na(cluster)) %>%
  #filter(region == "CB") %>%
  ggplot(aes(x = EV1, y = EV2, fill = Group))+
  geom_point(aes(fill = Group), colour = "black", size=5, pch = 21)+
  xlab(paste0("PC1 ", "(",var_prop[1], "%)"))+
  ylab(paste0("PC2 ", "(",var_prop[2], "%)"))+
  theme_base() +
  theme(legend.background = element_blank(),
        text = element_text(size = 14),
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        legend.justification = c(1, 1), 
        legend.position = c(1, 1),
        legend.position = c(1, 1),
        #axis.title.y = element_text(margin=margin(0,6,0,0), face = "bold"),
        #axis.title.x = element_text(margin=margin(12,0,0,0), face = "bold"),
        plot.background = element_rect(color = NA))+
  scale_fill_manual(values = c("#77AB43", "#008FD5", "#FFFFFF"))

# box plot of sizes

# quick fix for cluster labels
cluster_names <- c("Mainland", "Bras d'Or", "White")

sl_plot <- pca_df %>%
  mutate(Sex = sex) %>%
  mutate(Sex = ifelse(Sex == "M", "Male", "Female")) %>%
  mutate(Group = cluster_names[match(cluster, cluster_short)]) %>%
  #filter(!(id %in% c("SR20", "LN29", "SR15", "2012_SF16"))) %>%
  filter(!is.na(cluster)) %>%
  #filter(region == "CB") %>%
  ggplot(aes(x = Group, y = sl, fill = Group))+
  geom_jitter(aes(fill = Group), colour = "black", size=5, pch = 21 , alpha = 0.75)+
  stat_summary(fun.y=median, fun.ymin=median, fun.ymax=median, 
               geom="crossbar", width=0.9, color = "black", fatten = 6)+
  ylab("Standard length (cm)")+
  xlab("Genotypic cluster")+
  theme_base() +
  theme(legend.key = element_blank(),
        text = element_text(size = 14),
        legend.position = "none",
        #axis.title.y = element_text(margin=margin(0,12,0,0), face = "bold"),
        #axis.title.x = element_text(margin=margin(12,0,0,0), face = "bold"),
        #axis.text.y = element_text(angle = 90, hjust = 0.5, margin=margin(0,8,0,0)),
        plot.margin=unit(c(0,0,0,0),"mm"),
        plot.background = element_rect(color = NA))+
  scale_fill_manual(values = c("#77AB43", "#008FD5", "#FFFFFF"))

pdf(file = "figures/FigureS1.pdf", width = 8.5, height = 8.5, useDingbats=FALSE)
       
# create an apporpriate viewport.  Modify the dimensions and coordinates as needed
vp.top_right <- viewport(height=unit(.5, "npc"), width=unit(0.5, "npc"), 
                           just=c("left","top"), y=1, x=0.5)

vp.top_left<- viewport(height=unit(.5, "npc"), width=unit(0.5, "npc"), 
                         just=c("left","top"), y=1, x=0)
            
# plot your base graphics 

#layout(matrix(c(1,2,3,4), 2, 2, byrow = TRUE), widths=c(1,1), heights=c(1,1))
margins_maps <- c(4,4,1,1)
par(mfrow=c(2,2), mar = margins_maps)
frame()
frame()
                 
# plot the ggplot using the print command


# whole provence 
plot(1, xlim=c(-66.5,-59.75), ylim=c(43.25,47.25), 
     ylab="Longitude", xlab="Latitude", axes = FALSE, mar= margins_maps, cex.lab = 1.2)
map("worldHires","Canada", xlim=c(-66.5,-59.75), ylim=c(43.25,47.25), col="gray90", fill=TRUE,
    resolution = 0, mar = margins_maps, ylab = "Latitude", ylab = "Longitude", add = TRUE)
map.axes(cex.axis = 1.2)

lapply(1:length(pop.coords[[1]]), plot_pies_cluster, pop.coords = pop.coords, radius = 0.1) %>%  invisible()
lapply(1:length(pop.coords[[1]]), plot_labels,  pop.dat =  pop.dat) %>%  invisible()

rect(xleft = -62.125, ybottom = 45.25, xright = -60.5, ytop = 46.2)
map.scale(x = -62, y= 44, ratio = FALSE, cex = 0.7 )


# straight of canso region

plot(1, xlim=c(-62.125,-60.5), ylim=c(45.25,46.2), 
     ylab="Longitude", xlab="Latitude", axes = FALSE, mar= margins_maps, cex.lab = 1.2)
map("worldHires","Canada", xlim=c(-62.125,-60.5), ylim=c(45.25,46.2), col="gray90", fill=TRUE,
    resolution = 0, margins= margins_maps, add = TRUE)
#axis(side = 1)
map.axes(cex.axis = 1.2)

lapply(1:length(pop.coords[[1]]), plot_pies_cluster, pop.coords = pop.coords, radius = 0.025) %>%  invisible()
lapply(1:length(pop.coords[[1]]), plot_labels,  pop.dat =  pop.dat, zoom = TRUE) %>%  invisible()

map.scale(x = -60.85, y= 45.4, ratio = FALSE, cex = 0.7 )

vps <- baseViewports()
print(sl_plot, vp=vp.top_right)
print(pca_clust, vp=vp.top_left)


dev.off()


