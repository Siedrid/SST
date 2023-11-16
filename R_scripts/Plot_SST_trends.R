# Plot SST Trend Maps

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

devtools::install_github('oswaldosantos/ggsn')

# Read MK Datasets
setwd("E:/TIMELINE_SST/")
it_mosaics <- "OUT/Mosaics/"
plt_path = "Maps/"
#shp_bs <- "GIS/coast_dist/Baltic_Sea_coast.shp"
shp_path <- "GIS/Europe/Europe.gpkg"
africa_path <- "GIS/Africa/afr_g2014_2013_0.shp"
study_site <- "GIS/sst_analysis_polygons/study_area.shp"
poly_path <- "GIS/sst_analysis_polygons/intersting_sst_analysis.shp"
site_path <- "GIS/Sites/"
ocean_path <- "GIS/World_Seas_IHO_v3/"

# Read shapefiles
study_shp <- sf::st_read(study_site) # Polygon displaying boarders of study area
sf_object <- sf::st_read(shp_path) # Europe 
africa_shp <- sf::st_read(africa_path) # Africa
poly_shp <- sf::st_read(poly_path) # grid of tiles
ocean_shp <- sf::st_read(ocean_path)

europe_crs <- st_crs("+init=epsg:3035")

# Transform Shapefiles
baltic <- st_transform(sf_object, europe_crs)
study_shp <- st_transform(study_shp, europe_crs)
africa_shp <- st_transform(africa_shp, europe_crs)
ocean_shp <- st_transform(ocean_shp, europe_crs)
poly_shp <- st_transform(poly_shp, europe_crs)

# load and transform shapefiles
load_shp <- function(shp_path){
  europe_crs <- st_crs("+init=epsg:3035")
  shp <- sf::st_read(shp_path) %>% st_transform(europe_crs)
  return(shp)
}
# save grid in the end
save_grid <- function(grob_lst, map_name){
  g <- marrangeGrob(grobs = grob_lst, ncol =2, nrow=3)
  ggsave(paste0(plt_path, map_name), g, width=20, height = 30, units = "cm")
  dev.off()
}

# create SpatRaster
stack_lst <- function(short){
  r_list <- list()
  for (m in 1:12){
    m_str <- sprintf("%02d", m)
    mosaic_name <- paste0(it_mosaics, short, '/', m_str, "_merged_mosaic_mk_", short, ".nc")
    r_list[[m]] <- terra::rast(mosaic_name, lyrs = 'slope')
  }
  r_c <- terra::rast(r_list)
  names(r_c) <- month.abb
  return(r_c)
}


mask_IHO <- function(ocean_name, sp_raster, bb){
  #poly_ids <- st_crop(poly_shp, bb)
  IHO <- crop(sp_raster, bb)
  IHO <- mask(IHO, roi)
  return(IHO)
}

# background shapefile for remaining ocean not in study area
crop_ocean <- function(ocean_name){
  oc <- ocean_shp[!ocean_shp$NAME %in% ocean_name,] %>% st_crop(bb)
  oc <- oc[oc$NAME != 'South Pacific Ocean',]
  return(oc)
}

ocean_name = 'Balearic (Iberian Sea)'
ocean_name = c('Skagerrak', 'Kattegat', 'Baltic Sea')
ocean_name = c('Sea of Marmara', 'Aegean Sea') # Greece
ocean_name = 'Black Sea'
ocean_name = 'North Sea'

roi <- ocean_shp %>% filter(NAME %in% ocean_name)
bb <- st_bbox(roi)
short = 'baltic'

r_c <- stack_lst(short)
masked_rast <- mask_IHO(ocean_name, r_c,bb)
oc <- crop_ocean(ocean_name)

g <- ggplot()+
  geom_spatraster(data = masked_rast)+
  facet_wrap(~lyr, ncol = 3)+
  scale_fill_gradientn(name = 'K/year', colours = rev(brewer.pal(9, 'Spectral')), na.value = 'white', limits = c(-0.1, 0.1))+
  new_scale_fill()+
  #geom_sf(data = medi, aes(fill = 'grey'))+
  geom_sf(data = baltic, aes(fill = 'light grey')) +
  geom_sf(data = oc, aes(fill = 'grey'))+
  #geom_sf(data = roi, fill = NA)+
  scale_fill_manual(name = "", values = c("light grey", "dark grey"), label = c("Sea outside study area", "Land"))+
  coord_sf(xlim = c(bb[1],bb[3]), ylim=c(bb[2],bb[4]), expand = F)+
  xlab("Longitude")+
  ylab("Latitude")+
  ggtitle("SST anomaly trends")+
  scale_y_continuous(breaks = seq(34,65, by = 1))+
  scale_x_continuous(breaks = seq(22,40, by = 2))+ # only for Greece
  theme(panel.grid.major = element_line(color = gray(.5), linetype = "dashed", 
                                        linewidth = 0.5), 
        panel.ontop = T,
        panel.background = element_rect(fill = NA),
        strip.background = element_rect(fill = "white"),
        legend.position = "bottom", legend.key.width = unit(0.5, "inches"))

g
ggsave(paste0(plt_path, 'NorthSea-v1.png'), g, width = 20, height = 30, units = "cm")


# Read MK Mosaics with ocean shape
sst_trends <- function(short, ocean_name){
  m_lst <- c(1:12)
  
  # get lst of tile IDs covering study site .gpkg
  roi <- ocean_shp[ocean_shp$NAME == ocean_name,]
  medi <- ocean_shp[ocean_shp$NAME == "Mediterranean Sea - Western Basin",]
  
  # get extents of study site
  bb <- st_bbox(roi)
  
  for (m in m_lst){
    m = 1
    m_str <- sprintf("%02d", m)
    v_name <- paste0("gg_",m_str)
    mosaic_name <- paste0(it_mosaics, short, '/', m_str, "_merged_mosaic_mk_", short, ".nc")
    spain_mosaic <- paste0(it_mosaics, 'malle/', m_str, "_merged_mosaic_mk_malle.nc") 
    # for site Spain, two tiles are missing in first mosaic
    stars_object <- mask_IHO(ocean_name = ocean_name, mosaic_name = mosaic_name) %>% st_as_stars()
    #stars_object <- mask(terra::rast(mosaic_name, lyrs = 'slope'),roi) %>% st_as_stars()
    #malle <- mask(terra::rast(spain_mosaic, lyrs = 'slope'),roi) %>% st_as_stars()
    malle <- mask_IHO(ocean_name = ocean_name, mosaic_name = spain_mosaic) %>% st_as_stars()
    
    
    ggplot()+
      #scale_fill_manual(name = "Land cover", values = c("))+
      geom_stars(data = stars_object)+
      geom_stars(data = malle) +
      scale_fill_gradientn(name = 'K/year', colours = rev(brewer.pal(9, 'Spectral')), na.value = 'white', limits = c(-0.1, 0.1))+
      new_scale_fill()+
      geom_sf(data = medi, aes(fill = 'grey'))+
      #geom_sf(data = study_shp$geometry, fill = NA)+
      geom_sf(data = baltic, aes(fill = 'light grey')) +
      #geom_sf(data = roi, fill = NA)+
      scale_fill_manual(name = "Land cover", values = c("light grey", "dark grey"), label = c("Sea outside study area", "Land"))+
      coord_sf(xlim = c(bb[1],bb[3]), ylim=c(bb[2],bb[4]), expand = F)+
      xlab("Longitude")+
      ylab("Latitude")+
      ggtitle(month.name[m])+
      theme(panel.grid.major = element_line(color = gray(.5), linetype = "dashed", 
                                            linewidth = 0.5), 
            panel.ontop = T,
            panel.background = element_rect(fill = NA),
            legend.position = "bottom", legend.key.width = unit(0.1, "inches"), legend.direction = "vertical")
      #scale_fill_(name="Land", values = c("black", "blue"))
    
    
    assign(v_name, plt)
    
  }
  res <- list(gg_01, gg_02, gg_03, gg_04, gg_05, gg_06, gg_07, gg_08, gg_09, gg_10, gg_11, gg_12)
  #save_grid(grob_lst = res, strsplit(shape_name, ".", fixed = T)[[1]][1])
  return(res)
}

res <- sst_trends('spain', ocean_name)

g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

plt2 <- ggplot()+
  geom_sf(data = study_shp$geometry, aes(fill = "grey"))+
  geom_sf(data = baltic, aes(fill = "light grey"), show.legend = 'polygon') +  
  scale_fill_manual(name = "Land cover", labels = c("land", "sea"), values = c("red", "blue"))

legend_2 <- g_legend(plt2)
legend_1 <- g_legend(plt)
#coord_sf(xlim = c(bb[1],bb[3]), ylim=c(bb[2],bb[4]), expand = F)
library(ggpubr)
library(cowplot)
grid1 <- ggarrange(res[[1]] + theme(legend.position = "none"), res[[2]] + theme(legend.position = "none"))
legends <- ggarrange(legend_1, legend_2, heights = 0.2)
ggarrange(grid1, legends, nrow = 2)

# Read MK Mosaics
sst_trends <- function(short, shape_name){
  m_lst <- c(1:12)
  m = 1
  # get lst of tile IDs covering study site .gpkg
  site_shp <- sf::st_read(paste0(site_path, shape_name))
  poly_lst <- as.character(site_shp$id)
  
  # get extents of study site
  study_site <- poly_shp[poly_shp$id %in% poly_lst, ]
  site_tr <- st_transform(study_site, europe_crs)
  bb <- st_bbox(site_tr)
  
  for (m in m_lst){
    m = 1
    m_str <- sprintf("%02d", m)
    v_name <- paste0("gg_",m_str)
    mosaic_name <- paste0(it_mosaics, short, '/', m_str, "_merged_mosaic_mk_", short, ".nc")
    spain_mosaic <- paste0(it_mosaics, 'malle/', m_str, "_merged_mosaic_mk_malle.nc") 
    # for site Spain, two tiles are missing in first mosaic
    stars_object <- terra::rast(mosaic_name, lyrs = 'slope') %>% st_as_stars()
    malle <- terra::rast(spain_mosaic, lyrs = 'slope') %>% st_as_stars()
    
    plt <- ggplot()+  
      geom_stars(data = stars_object)+
      geom_stars(data = malle) +
      scale_fill_gradientn(name = 'K/year', colours = rev(brewer.pal(9, 'RdBu')), na.value = NA, limits = c(-0.1, 0.1))+
      geom_sf(data = study_shp, fill = NA, color = "grey")+
      geom_sf(data = baltic) + 
      coord_sf(xlim = c(bb[1],bb[3]), ylim=c(bb[2],bb[4]), expand = F)+
      xlab("Longitude")+
      ylab("Latitude")+
      ggtitle(month.name[m])+
      theme(panel.grid.major = element_line(color = gray(.5), linetype = "dashed", 
                                            linewidth = 0.5), panel.background = element_rect(fill = "white"),
            legend.position = "right", legend.key.width = unit(0.1, "inches"))
    
    
    assign(v_name, plt)
    
  }
  res <- list(gg_01, gg_02, gg_03, gg_04, gg_05, gg_06, gg_07, gg_08, gg_09, gg_10, gg_11, gg_12)
  save_grid(grob_lst = res, strsplit(shape_name, ".", fixed = T)[[1]][1])
  return(res)
}
short = "spain"
shape_name = "Malle.gpkg"

res <- sst_trends("spain", 'Malle.gpkg')
res_baltic <- sst_trends("baltic", "Denmark.gpkg")
res_italy <- sst_trends("italy", "Italy.gpkg")
res_greece <- sst_trends("eastern_europe", "Greece.gpkg")
res_black_sea <- sst_trends("eastern_europe", "Black_Sea.gpkg") 
res[1]

