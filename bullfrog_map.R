###Map ponds
#load packages
library(ggplot2)
library(maps)
library(sf)
library("rnaturalearth")
library("rnaturalearthdata")
library('usmap')
library(ggpubr)

##Map of the world! with isolates
#helpful tutorial: https://r-spatial.org/r/2018/10/25/ggplot2-sf-2.html
#read GPS
gps<-read.delim("~/Documents/GitHub/Bullfrog-nutrient-enrichment/gps.txt")

#load in world map
world <- ne_countries(scale = "medium", returnclass = "sf")

#get state borders
states_of_interest <- c("MA","RI", "CT", 'NH', "VT")
states <- us_map(include = states_of_interest)

#plot mass
maplot<-ggplot(data=world)+
  geom_sf()+
  geom_point(data=gps, aes(y=Lat, x=Log), color='black', size=4)+
  scale_color_manual(values=c('blue', 'orange'))+
  ylab('Latitude')+
  theme_bw()+
  ylab("Longitude")+
  borders("state") +
coord_sf(ylim = c(40, 45), xlim = c(-70,-78), expand = FALSE)

#plot just the sites
siteplot<-ggplot(data=world)+
  geom_sf()+
  geom_point(data=gps, aes(y=Lat, x=Log, color=Type), size=4)+
  scale_color_manual(values=c('orange', 'blue'))+
  borders("state") +
  coord_sf(ylim = c(41.98, 42.04), xlim = c(-71.1, -71.18), expand = FALSE)

ggarrange(maplot, siteplot)


