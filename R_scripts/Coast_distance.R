library(raster)
library(terra)
library(sf)
library(RColorBrewer)
library(tidyr)
library(ggplot2)
library(dplyr)

get_short <- function(ocean_name){
  
  # returns short from IHO ocean database, used also to name the Mosaic .nc files in Mosaic.py
  lst <- list.files(tile_path)
  short <- gsub(" ", "_", ocean_name)
  short <- gsub(".csv", "", short)
  short <- gsub("[()]", "", short)
  
  return(short)
}

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
# load and transform shapefiles to epsg 3035
load_shp <- function(shp_path){
  shp <- sf::st_read(shp_path) %>% st_transform(europe_crs, crs = 3035)
  return(shp)
}

# create SpatRaster from MK datasets, layer = slope
stack_lst <- function(short){
  # Read Statistic of valid trends observations
  trend_stats <- read.csv(paste0("E:/Publications/SST_analysis/Stats/","stats_trend_",short,".csv"))
  threshold<-30
  r_list <- list()
  month_list <- c()
  for (m in 1:12){
    perc_valid<-as.numeric(trend_stats[trend_stats$short == short & trend_stats$month==m,]['percent'])
    if (perc_valid > threshold){
      m_str <- sprintf("%02d", m)
      mosaic_name <- paste0(mosaic_path, short, '/', m_str, "_merged_mosaic_mk_", short, ".nc")
      r_list[[m]] <- terra::rast(mosaic_name, lyrs = 'slope')
      month_list<-c(month_list,m)
    }
  }
  r_c <- terra::rast(r_list)
  names(r_c) <- month.abb[month_list]
  return(r_c)
}

mask_IHO <- function(sp_raster, bb, roi){
  IHO <- crop(sp_raster, bb)
  IHO <- mask(IHO, roi)
  return(IHO)
}

coast_distance_slope <- function(m_list){
  
  monthly_df <- data.frame(matrix(ncol = length(m_list), nrow = 50))
  colnames(monthly_df) <- month.abb

  europe_cropped <- sf::st_crop(europe, bb)
  
  #for (m in 1:length(m_list)){
  for (m in names(masked_rast)){
    
    print(paste(m, "beeing processed"))
    mean_slope <- rep(NaN, 50)
    
    for (i in c(1:50)){
      print(paste(i, "km buffer beeing processed"))
      buf_2 <- st_buffer(europe_cropped, 1000 * i)
      buf_1 <- st_buffer(europe_cropped, 1000 * (i-1))
    
      dif <- st_difference(buf_2, buf_1)
    
      masked_slope <- mask(masked_rast[[m]], dif)
      mean_slope[i] <- mean(masked_slope[,,1], na.rm = T)
    }
    monthly_df[,m] <- mean_slope
  }
  return(monthly_df)
}

write_coast_dist <- function(monthly_df){
  outfile <- paste0(path_out, short, '_coast_dist.csv')
  write.csv(monthly_df, outfile)
}

# Plots ----

load_coast_dist <- function(study_areas){
  short <- study_areas[[1]][1] %>% get_short()
  df_all_sites <- read.csv(paste0(path_out, short, "_coast_dist.csv"))
  df_all_sites$site <- rep(short, nrow(df_all_sites))
  df_all_sites$dist <- c(1:50)
  
  for (i in 2:length(study_areas)){
    s <- study_areas[[i]][1] %>% get_short()
    df <-read.csv(paste0(path_out, s, "_coast_dist.csv"))
    df$site <- rep(s, nrow(df))
    df$dist <- c(1:50)
    df_all_sites <- rbind(df, df_all_sites)
  }
  
  # Reshape the data using gather
  df_all_sites2 <- df_all_sites[,2:ncol(df_all_sites)]
  df_long <- gather(df_all_sites2, key = "Variable", value = "Value", -dist, -site)
  df_long$Month <- factor(df_long$Variable, levels = unique(df_long$Variable))
  
  return(df_long)
}

# Create a line plot for each column using facet_wrap
plot_coast_dist <- function(df_long){
  df_long$Value<-df_long$Value*10
  g <- ggplot(df_long, aes(x = dist, y = Value, color = Month)) +
  geom_line() +
  facet_wrap(~ site, scales = "free_y") +
  labs(title ="Senn's Slope vs. Coast Distance")+
  xlab("Distance [km]")+
  ylab("Slope [K/decade]")+
  guides(fill=guide_legend(title="Month"))+
  ylim(0, 1)
  ggsave(paste0(plt_path, 'coast_dist.png'))
}

# Main Workflow ----

# setwd("E:/TIMELINE_SST/")
# shp_path <- "GIS/Europe/Europe.gpkg"
# ocean_path <- "GIS/World_Seas_IHO_v3/"
# 
# plt_path = "OUT/Plots/Coast_distance/"
# path_out = "GIS/coast_dist/" # directory to write output csv files to
# mosaic_path <- "OUT/Mosaics/" #path_out von Mosaic.py mit Unterordnern "Skagerrak", Adriatic_Sea", etc.
# tile_path <- "E:/TIMELINE_SST/Tile_Lists/"

shp_path <- "E:/Publications/SST_analysis/GIS Projekt/Europe.gpkg"
ocean_path <- "E:/Conferences_Presentations/Strukturkommision_2022/Folien/World_Seas_IHO_v3/World_Seas_IHO_v3/"

plt_path = "E:/Publications/SST_analysis/Figures/Figure11/new_test/"
path_out = "E:/Publications/SST_analysis/Figures/Figure11/new_test/"
#mosaic_path <- "E:/Publications/SST_analysis/Mosaics/New/"
mosaic_path <- "E:/Publications/SST_analysis/Mosaics/cropped_mk/"
tile_path <- "E:/Publications/SST_analysis/to_process/"

A = c('Skagerrak', 'Kattegat', 'Baltic Sea', 'North Sea')
B = 'Adriatic Sea'
D = c('Aegean Sea', 'Sea of Marmara')
E = 'Balearic (Iberian Sea)'

study_areas <- list(A,B,D,E)

europe <- load_shp(shp_path)
ocean_shp <- load_shp(ocean_path)
m_list <- month.abb

for (i in 1:length(study_areas)){
  ocean_name <- study_areas[[i]]
  
  short <- get_short(ocean_name[1])
  
  tile_df <- read_tile_lists(tile_path, ocean_name[1])
  roi <- ocean_shp %>% filter(NAME %in% ocean_name)
  bb <- get_bb_list(tile_df)
  short <- get_short(ocean_name[1])
  
  #r_c <- stack_lst(short)
  #masked_rast <- mask_IHO(r_c, bb, roi)
  masked_rast <- stack_lst(short)
  
  monthly_df <- coast_distance_slope(m_list)
  write_coast_dist(monthly_df)
  
}

df_long <- load_coast_dist(study_areas)
plot_coast_dist(df_long)