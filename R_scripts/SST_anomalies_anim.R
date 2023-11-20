## Animations of SST anomalies across European Boarders

library(ggplot2)
library(tidyverse)
library(maps)
library(scales)
library(stringr)
library(terra)
library(stars)
library(RColorBrewer)
library(ggspatial)
library(gridExtra)
library(ggnewscale)
library(tidyterra)
#library(ggsn)
library(dplyr)
library(animation)

setwd("E:/TIMELINE_SST/")
it_mosaics <- "OUT/Mosaics/"
europe_path <- "GIS/Europe/Europe.gpkg"
africa_path <- "GIS/Africa/afr_g2014_2013_0.shp"

# load and transform shapefiles
load_shp <- function(shp_path){
  europe_crs <- st_crs("+init=epsg:3035")
  shp <- sf::st_read(shp_path) %>% st_transform(europe_crs)
  return(shp)
}
euro <- load_shp(europe_path)
afr <- load_shp(africa_path)

short = 'eastern_europe'
month = '06'

short_lst = c('eastern_europe', 'spain', 'baltic', 'french_atlantic', 'GB', 'italy', 'malle', 'northern_africa')

m_lst <- c(1:12)

animate_sst_anoms <- function(m){
  m_str <- sprintf("%02d", m)
  
  for (i in 1:length(short_lst)){
    sst <- terra::rast(paste0(it_mosaics, short_lst[i], '/', m_str, '_merged_mosaic_dif_2022_', short_lst[i], '.nc'), lyrs = 'sst_dif_max')
    assign(paste0('sst', i), sst)
  } 
  plt <- ggplot()+
    geom_spatraster(data = sst1)+
    geom_spatraster(data = sst2)+
    geom_spatraster(data = sst3)+
    geom_spatraster(data = sst4)+
    geom_spatraster(data = sst5)+
    geom_spatraster(data = sst6)+
    geom_spatraster(data = sst7)+
    geom_spatraster(data = sst8)+
    
    scale_fill_gradientn(name = 'K', colours = rev(brewer.pal(9, 'Spectral')), na.value = NA, limits = c(-3, 3))+
    geom_sf(data = euro)+
    geom_sf(data = afr)+
    coord_sf(xlim = c(2521500,6754500), ylim=c(900500,4516500), expand = F)+
    ggtitle('SST anomalies with respect to the period 1990-2022', subtitle = paste0('2022-', m_str))+
    theme(panel.grid.major = element_line(color = gray(.5), linetype = "dashed", 
                                          linewidth = 0.5), 
          panel.ontop = T,
          panel.background = element_rect(fill = NA),
          strip.background = element_rect(fill = "white"),
          legend.position = "bottom", legend.key.width = unit(0.5, "inches"))
  
}

my_plots <- lapply(m_lst, animate_sst_anoms)

saveGIF({
  for (i in 1:12) plot(my_plots[[i]])},
  movie.name = "SST_2022.gif"
)
