library(raster)
library(terra)
library(sf)
library(RColorBrewer)
library(tidyr)
library(ggplot2)
# edited 
setwd("E:/TIMELINE_SST/")
shp_path <- "GIS/Europe/Europe.gpkg"
plt_path = "OUT/Plots/Coast_distance/"
mosaic_path <- "OUT/Mosaics/" #path_out von Mosaic.py mit Unterordnern "Skagerrak", Adriatic_Sea", etc.
tile_path <- "E:/TIMELINE_SST/Tile_Lists/"

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

# create SpatRaster
stack_lst <- function(short){
  r_list <- list()
  for (m in 1:12){
    m_str <- sprintf("%02d", m)
    mosaic_name <- paste0(mosaic_path, short, '/', m_str, "_merged_mosaic_mk_", short, ".nc")
    r_list[[m]] <- terra::rast(mosaic_name, lyrs = 'slope')
  }
  r_c <- terra::rast(r_list)
  names(r_c) <- month.abb
  return(r_c)
}


mask_IHO <- function(sp_raster, bb, roi){
  #poly_ids <- st_crop(poly_shp, bb)
  IHO <- crop(sp_raster, bb)
  IHO <- mask(IHO, roi)
  return(IHO)
}

europe <- load_shp(shp_path)
ocean_shp <- load_shp(ocean_path)

ocean_name <- 'Adriatic Sea'
short <- get_short(ocean_name)

tile_df <- read_tile_lists(tile_path, ocean_name[1])
roi <- ocean_shp %>% filter(NAME %in% ocean_name)
bb <- get_bb_list(tile_df)
short <- get_short(ocean_name[1])

r_c <- stack_lst(short)
masked_rast <- mask_IHO(r_c, bb, roi)

coast_distance_slope <- function(m_list, shape, short){
  
  path2mosaic <- paste0(mosaic_path, short, "/")
  monthly_df <- data.frame(matrix(ncol = length(m_list), nrow = 50))
  colnames(monthly_df) <- month.abb

  europe <- load_shp(shp_path)
  europe_cropped <- sf::st_crop(europe, bb)
  
  for (m in 1:length(m_list)){
    
    print(paste(m_list[m], "beeing processed"))
    mean_slope <- rep(NaN, 50)
    
    for (i in c(1:50)){
      buf_2 <- st_buffer(europe_cropped, 1000 * i)
      buf_1 <- st_buffer(europe_cropped, 1000 * (i-1))
    
      dif <- st_difference(buf_2, buf_1)
    
      masked_slope <- mask(masked_rast[[i]], dif)
      mean_slope[i] <- mean(masked_slope[,,1], na.rm = T)
    }
    monthly_df[,m] <- mean_slope
  }
  
  return(monthly_df)
}

write_coast_dist <- function(monthly_df){
  outfile <- paste0(path_out, short, 'coast_dist.csv')
  write.csv(monthly_df, outfile)
}

m_list <- c("01", "02", "03", "06", "07", "09", "12")

# Black Sea
df <- distance_slope(m_list, shape, "eastern_europe")
write.csv(df, "GIS/coast_dist/Black_sea_coast_dist-v3.csv")

# Baltic
baltic_df <- distance_slope(m_list, "Denmark.gpkg", "baltic")
write.csv(baltic_df, "GIS/coast_dist/Denmark_coast_dist-v3.csv")

# Italy
it_mosaics <- "OUT/Mosaics/italy/"
italy_df <- distance_slope(m_list, "Italy.gpkg", "italy")
write.csv(italy_df, "GIS/coast_dist/Italy_coast_dist-v3.csv")

# Greece
greece_df <- distance_slope(m_list, "Greece.gpkg", "eastern_europe")
write.csv(greece_df, "GIS/coast_dist/Greece_coast_dist-v3.csv")

# Plots ----


plot_coast_distance <- function(coast_df, site){
  coast_dist <- seq(1, 50, by = 1)
  #colPal = viridis_pal(option = "D")(length(greece_df))
  colPal = brewer.pal(n=length(coast_df), name = "Dark2")
  #png(paste0(plt_path, site, '-v4.png'),  width = 6, height = 4, units = "in", res = 1200)
  par(mar = c(4, 4.1, 2, 2.1))
  plot(coast_dist, coast_df[,1], xlab = "Coast Distance [km]", ylab = "Slope [K/year]", col = colPal[1], ylim = c(0.03,0.085), pch = 16)
  lines(coast_dist, coast_df[,1], lwd = 2, col = colPal[1])
  
  for (i in 2:length(coast_df)){
    points(coast_dist, coast_df[,i], col = colPal[i], pch =16)
    lines(coast_dist, coast_df[,i], col = colPal[i], lwd =2)
  }
  grid()
  legend_labs <- colnames(coast_df)
  
  legend("topright", legend_labs, col = colPal, lwd = 2)
  title(paste("Senn's Slope vs. Coast Distance for", site))
  #dev.off()
}

# load Dataframes and concatenate
site = 'Black_Sea'
df_all_sites <- read.csv(paste0("GIS/coast_dist/", site, "_coast_dist-v3.csv"))
df_all_sites$site <- rep(site, nrow(df_all_sites))
df_all_sites$dist <- c(1:50)

site_lst <- c("Denmark", "Greece", "Italy")

for (s in site_lst){
  df <-read.csv(paste0("GIS/coast_dist/", s, "_coast_dist-v3.csv"))
  df$site <- rep(s, nrow(df))
  df$dist <- c(1:50)
  #df <- df[2:10]
  df_all_sites <- rbind(df, df_all_sites)
}

# Reshape the data using gather
#df$dist <- c(1:50)
df_all_sites2 <- df_all_sites[,2:10]
colnames(df_all_sites2) <- c("Jan", "Feb", "Mar", "Jun", "Jul", "Sep", "Dec", "site", "dist")
df_long <- gather(df_all_sites2, key = "Variable", value = "Value", -dist, -site)
df_long$Month <- factor(df_long$Variable, levels = c("Jan", "Feb", "Mar", "Jun", "Jul", "Sep", "Dec"))


# Create a line plot for each column using facet_wrap
x11()
png(paste0(plt_path, 'allsites_coast_dist-v1.png'),  width = 6, height = 4, units = "in", res = 1200)
ggplot(df_long, aes(x = dist, y = Value, color = Month)) +
  geom_line() +
  facet_wrap(~ site, scales = "free_y") +
  labs(title ="Senn's Slope vs. Coast Distance")+
  xlab("Distance [km]")+
  ylab("Slope [K/year]")+
  guides(fill=guide_legend(title="Month"))
dev.off()
