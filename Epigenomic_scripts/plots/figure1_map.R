#!/usr/bin/Rscript
########################################################################
#############                MAP FIGURE               ##################
########################################################################



library("rnaturalearth")
library("rnaturalearthdata")
library(ggplot2)
world <- ne_countries(scale = "medium", returnclass = "sf")
class(world)
library("ggspatial")
ggplot(data = world) +
  theme(panel.grid.major = element_line(color = "gray97",  size = 0.3), 
        panel.background = element_rect(fill = "white"),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())+
  geom_sf(color = "grey50", fill = "grey69") +
  
  #annotation_scale(location = "bl", width_hint = 0.5) + #because you have unprojected data, not to use
  annotation_north_arrow(location = "tl", which_north = "true",
                         #pad_x = unit(0.25, "in"), pad_y = unit(5.8, "in"),
                         style = north_arrow_fancy_orienteering) +
  
  annotate(geom = "text", x = 5.8 , y = 44.5, label = "Montpellier",
           color = "#5E3C99", size = 5) +
  annotate(geom = "point", x = 3.8772300, y = 43.6109200,
           color = "#5E3C99", size = 3) +
  annotate(geom = "text", x = 22.6 , y = 53.1, label = "Warsaw",
           color = "#FDB863", size = 5) +
  annotate(geom = "point", x = 21.0122, y = 52.2297 ,
           color = "#FDB863", size = 3) +
  annotate(geom = "text", x = 4.2 , y = 40.7, label = "Barcelona",
           color = "#E66101", size =5) +
  annotate(geom = "point", x = 2.1734, y = 41.3851 ,
           color = "#E66101", size = 3) +
  
  coord_sf(xlim = c(-12, 35), ylim = c(35, 60), expand=FALSE) +
  xlab("Longitude") + ylab("Latitude")

