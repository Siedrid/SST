library(raster)
library(terra)
library(sf)
library(viridis)
library(RColorBrewer)
library(tidyr)
library(ggplot2)
# edited 
setwd("E:/TIMELINE_SST/")
shp_path <- "GIS/Europe/Europe.gpkg"
plt_path = "OUT/Plots/Coast_distance/"

sites_path <- "GIS/Sites/"

# Site names are Black_Sea, Denmark, Greece, Italy which stand for the Shapefiles (.gpkg) in sites_path
# short refer to the shorts used for the different Mosaics
# shorts = c(baltic, eastern_europe, italy)

distance_slope <- function(m_list, shape, short){
  # mosaic_file: "_merged_mosaic_mk_baltic.nc"
  it_mosaics <- paste0("OUT/Mosaics/", short, "/")
  monthly_df <- data.frame(matrix(ncol = length(m_list), nrow = 50))
  colnames(monthly_df) <- m_list #month.abb
  europe_crs <- st_crs("+init=epsg:3035")
  
  europe <- st_read(shp_path)
  site <- st_read(paste0(sites_path, shape))
  site_tr <- st_transform(site, europe_crs)
  
  italy_transformed <- st_transform(europe, europe_crs)
  europe_cropped <- sf::st_crop(italy_transformed, st_bbox(site_tr))
  
  for (m in 1:length(m_list)){
    
    print(paste(m_list[m], "beeing processed"))
    mosaic_name <- paste0(m_list[m], "_merged_mosaic_mk_", short, ".nc")
    
    it_mk <- rast(paste0(it_mosaics,mosaic_name))

    mean_slope <- rep(NaN, 50)
    
    for (i in c(1:50)){
      buf_2 <- st_buffer(europe_cropped, 1000 * i)
      buf_1 <- st_buffer(europe_cropped, 1000 * (i-1))
    
      dif <- st_difference(buf_2, buf_1)
    
      maske_slope <- mask(it_mk$slope, dif)
      mean_slope[i] <- mean(maske_slope[,,1], na.rm = T)
    }
    monthly_df[,m] <- mean_slope
  }
  
  return(monthly_df)
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
