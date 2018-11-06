

### chlamy map
library(tidyverse)
library(readxl)
library(ggplot2)
library(viridis)
library(ggalt)
library(dplyr)
library(rgdal)
library(ggmap)
library(tools)
library(marmap)
library(dplyr)

locations <- read_excel("data-general/chlamy-population-key.xlsx")

loc2 <- locations %>% 
	mutate(lat = ifelse(is.na(lat), 42.37580, lat)) %>% 
	mutate(long = ifelse(is.na(long), -72.519867, long)) %>% 
	dplyr::rename(latitude = lat) %>% 
	dplyr::rename(longitude = long)

library(colormap)
library(sf)
library(maps)
library(viridis)
library(tidyverse)

global_therm <- read_csv("data-raw/global_therm.csv") %>% 
	filter(!is.na(Mean), Biome == "terrestrial") 


world1 <- sf::st_as_sf(map('world', plot = FALSE, fill = TRUE, col = 1:10, wrap=c(-180,180)))

ic <- colormap(colormap = colormaps$viridis, nshades = 8, format = "hex",
			   alpha = 1, reverse = FALSE)

topt_data <- st_as_sf(topts2, coords = c("longitude", "latitude"), crs = 4326)
topt_data <-st_transform(x = topt_data, crs = "+proj=robin")

global_data <- st_as_sf(global_therm, coords = c("Lon", "Lat"), crs = 4326)

global_data <-st_transform(x = global_data, crs = "+proj=robin")

topt_data <- st_as_sf(loc2, coords = c("longitude", "latitude"), crs = 4326)
topt_data <-st_transform(x = topt_data, crs = "+proj=robin")

ggplot(global_data) +
	geom_sf(aes(color = Mean)) +
	scale_color_viridis(discrete = FALSE, option = "inferno") +
	geom_sf(data = world1, color = "black", fill = "darkgrey") +
	# geom_sf(data = topt_data, geom = "point", aes(color = sst_sd), size = 3)+
	geom_sf(data = topt_data, geom = "point", shape = 1, color = "black", size = 3)+
	theme_bw()+
	theme(panel.grid = element_blank(),
		  line = element_blank(),
		  rect = element_blank(),
		  text = element_text(size=14))
ggsave("figures/global_map.pdf", width = 6.5, height = 4)


ggplot(global_therm,aes(x=Lon,y=Lat,fill=Mean))+
	geom_tile()+  
	geom_hline(yintercept = c(-60,-30,0,30,60),color="black",size=0.1)+
	geom_vline(xintercept = seq(-180,180,by=60),color="black",size=0.1)+
	scale_alpha(trans="log")+
	scale_fill_gradientn(colors=c("midnightblue","dodgerblue","#4dac26","gold" ,"#d7191c"))+
	theme(axis.line        = element_blank(),
		  axis.text        = element_blank(),
		  axis.ticks       = element_blank(),
		  axis.title       = element_blank(),
		  panel.background = element_blank(),
		  panel.border     = element_blank(),
		  panel.margin = unit(0,"null"),
		  panel.grid.major = element_blank(),
		  panel.grid.minor = element_blank(),
		  plot.background  = element_blank(),
		  plot.margin = rep(unit(0,"null"),4))+
	geom_sf(data = topt_data, geom = "point", shape = 1, color = "black", size = 3)+
	geom_contour(data = filter(coastline,y>-60),aes(x=x,y=y,z=z,fill=NULL),breaks=0, colour="grey40")+
	# geom_point(data=filter(all_Thermal,!is.na(Latitude)),aes(x=Longitude,y=Latitude,fill=TOpt,pch=Dataset,size=T_breadth))+
	scale_shape_manual(values = c(21:22))



