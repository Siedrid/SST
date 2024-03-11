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
library(dplyr)

read_tile_lists <- function(tile_path, IHO_name){
  # get list of tile IDs covering IHO basin (from Drive)
  df <- read.csv(paste0(tile_path, IHO_name, '.csv'))
  return(df)
}

get_bb_list <- function(tile_df){
  # get bounding box from .csv files
  tile_bb <- st_bbox(c(xmin = min(tile_df$left), xmax = max(tile_df$right), 
                       ymax = max(tile_df$top), ymin = min(tile_df$bottom)), crs = 4326)
  tile_bb <- st_transform(st_as_sfc(tile_bb), crs = 3035) %>% st_bbox()
  return(tile_bb)
}

get_short <- function(ocean_name){
  
  # returns short from IHO ocean database, used also to name the Mosaic .nc files in Mosaic.py
  lst <- list.files(tile_path)
  short <- gsub(" ", "_", ocean_name)
  short <- gsub(".csv", "", short)
  short <- gsub("[()]", "", short)
  
  return(short)
}

# load and transform shapefiles to epsg 3035
load_shp <- function(shp_path){
  shp <- sf::st_read(shp_path) %>% st_transform(europe_crs, crs = 3035)
  return(shp)
}


# create SpatRaster
stack_lst <- function(short, mosaic, year, layer, stats){
  # by default the year 2022 is used for the anomaly rasters, and the respective file name
  # Mosaic stands for either "MK" or "Dif", MK ds with the layers "slope" and "p" respectively,
  # or Anomaly ds with layer "obs_count" and "sst_dif_max"
  r_list <- list()
  #for (m in 1:12){
  i=1
  for (m in c(12,1,2,3,4,5,6,7,8,9,10,11)){
    m_str <- sprintf("%02d", m)
    if (mosaic == 'MK'){
      mosaic_name <- paste0(mosaic_path, short, '/', m_str, "_merged_mosaic_mk_", short, "_",stats,".nc")
      #mosaic_name <- paste0(mosaic_path, short, '/', m_str, "_merged_mosaic_mk_", short, ".nc")
      r_list[[i]] <- terra::rast(mosaic_name, lyrs = layer)
    } else if (mosaic == "Dif"){
      mosaic_name <- paste0(mosaic_path, short, '/', m_str, "_merged_mosaic_dif_",year,"_", short, "_",stats,".nc")
      r_list[[m]] <- terra::rast(mosaic_name, lyrs = layer)
    } else if (mosaic == "CCI"){
      mosaic_name <- paste0(path_cci_composites, m_str, "_merged_mosaic_dif_",year,"_", short, ".nc")
      r_list[[m]] <- terra::rast(mosaic_name, lyrs = layer)
    } else {
      print("Variable mosaic must be either MK or Dif")
    }
    i=i+1
  }
  r_c <- terra::rast(r_list)
  #names(r_c) <- month.abb
  names(r_c) <- c("Dec","Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov") 
  return(r_c)
}


mask_IHO <- function(sp_raster, bb, roi){
  IHO <- crop(sp_raster, bb)
  IHO <- mask(IHO, roi)
  return(IHO)
}

# background shapefile for remaining ocean not in study area
crop_ocean <- function(ocean_name, bb){
  oc <- ocean_shp[!ocean_shp$NAME %in% ocean_name,] %>% st_crop(bb)
  #oc <- ocean_shp %>% st_crop(bb)
  #oc <- ocean_shp[!ocean_shp$NAME %in% ocean_name,]
  oc <- oc[oc$NAME != 'South Pacific Ocean',]
  return(oc)
}

plot_sst_trend_maps <- function(masked_rast,stats){
  masked_rast<-masked_rast*10
  g <- ggplot()+
    geom_spatraster(data = masked_rast)+
    facet_wrap(~lyr, ncol = 3)+
    # Raster Legend:
    scale_fill_gradientn(name = 'K/decade', colours = rev(brewer.pal(9, 'Spectral')), na.value = 'white', limits = c(-1, 1))+
    new_scale_fill()+
    geom_sf(data = europe, aes(fill = 'light grey')) +
    geom_sf(data = oc, aes(fill = 'white'))+
    geom_sf(data = oc, aes(fill = 'grey'))+
    # Vector Layer Legend:
    scale_fill_manual(name = "", values = c("light grey", "dark grey", "white"), 
                      label = c("Sea outside study area", "Land", "Not significant trend"),
                      na.value = "white")+
    coord_sf(xlim = c(bb[1],bb[3]), ylim=c(bb[2],bb[4]), expand = F)+
    xlab("Longitude")+
    ylab("Latitude")+
    ggtitle("SST Anomaly Trends")+
    scale_y_continuous(breaks = seq(30,65, by = 2))+
    scale_x_continuous(breaks = seq(-15,40, by = 2))+ 
    # Grid theme:
    theme(panel.grid.major = element_line(color = gray(.5), linetype = "dashed", 
                                          linewidth = 0.5), 
          panel.ontop = T,
          panel.background = element_rect(fill = NA),
          strip.background = element_rect(fill = "white"),
          legend.position = "bottom", legend.key.width = unit(0.35, "inches"))+ # change scalebar length
    annotation_scale(location = 'bl')+
    annotation_north_arrow(location = 'bl', height = unit(0.75, "cm"), width = unit(0.75,"cm"),
                           pad_x = unit(0.1, "in"), pad_y = unit(0.3, "in"))
  
  ggsave(paste0(plt_path, short, '_anomaly_trends',"_",stats,'.png'), g, width = 20, height = 30, units = "cm")
}

plot_sst_anomaly_maps <- function(masked_rast,year,stats){
  
  g <- ggplot()+
    geom_spatraster(data = masked_rast)+
    facet_wrap(~lyr, ncol = 3)+
    # Raster Legend:
    scale_fill_gradientn(name = 'K', colours = rev(brewer.pal(9, 'Spectral')), na.value = 'white', limits = c(-3,3))+
    new_scale_fill()+
    geom_sf(data = europe, aes(fill = 'light grey')) +
    geom_sf(data = oc, aes(fill = 'white'))+
    geom_sf(data = oc, aes(fill = 'grey'))+
    # Vector Layer Legend:
    scale_fill_manual(name = "", values = c("light grey", "dark grey", "white"), 
                      label = c("Sea outside study area", "Land", "No Data"),
                      na.value = "white")+
    coord_sf(xlim = c(bb[1],bb[3]), ylim=c(bb[2],bb[4]), expand = F)+
    xlab("Longitude")+
    ylab("Latitude")+
    ggtitle(paste0("SST anomalies in ",year))+
    scale_y_continuous(breaks = seq(30,65, by = 2))+
    scale_x_continuous(breaks = seq(-15,40, by = 2))+ 
    # Grid theme:
    theme(panel.grid.major = element_line(color = gray(.5), linetype = "dashed", 
                                          linewidth = 0.5), 
          panel.ontop = T,
          panel.background = element_rect(fill = NA),
          strip.background = element_rect(fill = "white"),
          legend.position = "bottom", legend.key.width = unit(0.35, "inches"))+ #change scalebar length
    annotation_scale(location = 'bl')+
    annotation_north_arrow(location = 'bl', height = unit(0.75, "cm"), width = unit(0.75,"cm"), 
                           pad_x = unit(0.1, "in"), pad_y = unit(0.2, "in")) 
  
  ggsave(paste0(plt_path, short,'_',year, '_anomalies',"_",stats,'.png'), g, width = 20, height = 30, units = "cm")
}

plot_sst_anomaly_maps_cci <- function(masked_rast,year){
  
  g <- ggplot()+
    geom_spatraster(data = masked_rast)+
    facet_wrap(~lyr, ncol = 3)+
    # Raster Legend:
    scale_fill_gradientn(name = 'K', colours = rev(brewer.pal(9, 'Spectral')), na.value = 'white', limits = c(-2,2))+
    new_scale_fill()+
    geom_sf(data = europe, aes(fill = 'light grey')) +
    geom_sf(data = oc, aes(fill = 'white'))+
    geom_sf(data = oc, aes(fill = 'grey'))+
    # Vector Layer Legend:
    scale_fill_manual(name = "", values = c("light grey", "dark grey", "white"), 
                      label = c("Sea outside study area", "Land", "No Data"),
                      na.value = "white")+
    coord_sf(xlim = c(bb[1],bb[3]), ylim=c(bb[2],bb[4]), expand = F)+
    xlab("Longitude")+
    ylab("Latitude")+
    ggtitle(paste0("SST anomalies in ",year))+
    scale_y_continuous(breaks = seq(30,65, by = 2))+
    scale_x_continuous(breaks = seq(-15,40, by = 2))+ 
    # Grid theme:
    theme(panel.grid.major = element_line(color = gray(.5), linetype = "dashed", 
                                          linewidth = 0.5), 
          panel.ontop = T,
          panel.background = element_rect(fill = NA),
          strip.background = element_rect(fill = "white"),
          legend.position = "bottom", legend.key.width = unit(0.35, "inches"))+ #change scalebar length
    annotation_scale(location = 'bl')+
    annotation_north_arrow(location = 'bl', height = unit(0.75, "cm"), width = unit(0.75,"cm"), 
                           pad_x = unit(0.1, "in"), pad_y = unit(0.2, "in")) 
  
  ggsave(paste0(plt_path, short,'_',year, '_anomalies_cci.png'), g, width = 20, height = 30, units = "cm")
}



plot_p_value_maps <- function(masked_rast,stats){
  #masked_rast<-masked_rast*10
  
  g <- ggplot()+
    geom_spatraster(data = masked_rast)+
    facet_wrap(~lyr, ncol = 3)+
    # Raster Legend:
    scale_fill_gradientn(name = 'p-value', colours = rev(brewer.pal(9, 'Spectral')), na.value = 'red', limits = c(0, 0.4),oob=squish)+
    #scale_fill_gradientn(name = 'p-value', colours = rev(brewer.pal(9, 'Spectral')), limits = c(0, 0.3))+
    new_scale_fill()+
    geom_sf(data = europe, aes(fill = 'light grey')) +
    #geom_sf(data = oc, aes(fill = 'white'))+
    geom_sf(data = oc, aes(fill = 'grey'))+
    # Vector Layer Legend:
    scale_fill_manual(name = "", values = c("light grey", "dark grey"), 
                      label = c("Sea outside study area", "Land"),
                      na.value = "white")+
    coord_sf(xlim = c(bb[1],bb[3]), ylim=c(bb[2],bb[4]), expand = F)+
    xlab("Longitude")+
    ylab("Latitude")+
    ggtitle("SST anomaly trends")+
    scale_y_continuous(breaks = seq(30,65, by = 2))+
    scale_x_continuous(breaks = seq(-15,40, by = 2))+ 
    # Grid theme:
    theme(panel.grid.major = element_line(color = gray(.5), linetype = "dashed", 
                                          linewidth = 0.5), 
          panel.ontop = T,
          panel.background = element_rect(fill = NA),
          strip.background = element_rect(fill = "white"),
          legend.position = "bottom", legend.key.width = unit(0.35, "inches"))+ # change scalebar length
    annotation_scale(location = 'bl')+
    annotation_north_arrow(location = 'bl', height = unit(0.75, "cm"), width = unit(0.75,"cm"),
                           pad_x = unit(0.1, "in"), pad_y = unit(0.3, "in"))
  
  ggsave(paste0(plt_path, short, '_p_values',"_",stats,'.png'), g, width = 20, height = 30, units = "cm")
}

# Main Workflow-----

# setwd("E:/TIMELINE_SST/")
# mosaic_path <- "OUT/Mosaics/" #path_out von Mosaic.py mit Unterordnern "Skagerrak", Adriatic_Sea", etc.
# plt_path = "Maps/" # Pfad, wo Plots gespeichert werden
# 
# # Shape Paths
# shp_path <- "GIS/Europe/Europe.gpkg"
# africa_path <- "GIS/Africa/afr_g2014_2013_0.shp"
# study_site <- "GIS/sst_analysis_polygons/study_area.shp"
# poly_path <- "GIS/sst_analysis_polygons/intersting_sst_analysis.shp"
# ocean_path <- "GIS/World_Seas_IHO_v3/"

plt_path = "E:/Publications/SST_analysis/Figures_test/daytime_corr/"
#mosaic_path <- "E:/Publications/SST_analysis/Mosaics/max_daytime_cor/"
mosaic_path <- "E:/Publications/SST_analysis/Mosaics/Median/"
#mosaic_path <- "E:/Publications/SST_analysis/Mosaics/cropped_mk/"
path_cci_composites<-"E:/Publications/SST_analysis/CCI_composites/"

# Shape Paths
shp_path <- "E:/Publications/SST_analysis/GIS Projekt/Europe.gpkg"
africa_path <- "E:/Publications/SST_analysis/GIS Projekt/Africa/Africa/afr_g2014_2013_0.shp"
study_site <- "E:/Publications/SST_analysis/GIS Projekt/study_area/study_area.shp"
poly_path <- "E:/SST_Analysis/Shapes/intersting_sst_analysis.shp"
#ocean_path <- "E:/Conferences_Presentations/Strukturkommision_2022/Folien/World_Seas_IHO_v3/World_Seas_IHO_v3/"
ocean_path <- "E:/Publications/SST_analysis/Final_Study_Areas/World_Seas_Cropped/"
final_path <- "E:/Publications/SST_analysis/Final_Study_Areas/Final_areas/"


# path to tile lists
tile_path <- "E:/Publications/SST_analysis/to_process/"

# Transform Shapefiles
europe <- load_shp(shp_path)
africa_shp <- load_shp(africa_path)
ocean_shp <- load_shp(ocean_path)
poly_shp <- load_shp(poly_path)
final_area <- load_shp(final_path)

A = c('Skagerrak', 'Kattegat')
B = 'Adriatic Sea'
D = c('Aegean Sea', 'Sea of Marmara')
E = 'Balearic (Iberian Sea)'

#study_areas <- list(A,B,D,E)
study_areas <- list(B,E,D)

#year <- "1994"

years <- as.character(seq(from = 1995, to = 2022, by = 1))
#years <- list(1999)

stats<-'median'

for (year in years){

  for (i in 1:length(study_areas)){
    ocean_name <- study_areas[[i]]
    
    tile_df <- read_tile_lists(tile_path, ocean_name[1])
    #roi <- ocean_shp %>% filter(NAME %in% ocean_name)
    roi <-final_area
    bb <- get_bb_list(tile_df)
    xmin <-bb['xmin']
    oc <- crop_ocean(ocean_name,bb)
    
    short <- get_short(ocean_name[1])
    
    #r_c_mk <- stack_lst(short, 'MK', layer = 'slope',year, stats)
    #masked_rast_mk <- mask_IHO(r_c_mk, bb, roi)
    
    #r_c_p <- stack_lst(short, 'MK', layer = 'p')
    #masked_rast_p <- mask_IHO(r_c_p, bb, roi)
    
    r_c_dif <- stack_lst(short, 'Dif', year, layer = 'sst_dif_med', stats)
    masked_rast_dif <- mask_IHO(r_c_dif, bb, roi)
    
    #r_c_cci <- stack_lst(short, 'CCI', year, layer = 'sst_cci_anom')
    #masked_rast_cci <- mask_IHO(r_c_cci, bb, roi)
    
    
    #plot_sst_trend_maps(masked_rast_mk, stats)
    plot_sst_anomaly_maps(masked_rast_dif,year, stats)
    #plot_sst_anomaly_maps_cci(masked_rast_cci,year, stats)
    #plot_p_value_maps(masked_rast_p, stats)
  }
}