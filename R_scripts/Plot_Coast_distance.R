# Plot Coast distance

library(raster)
library(terra)
library(sf)
library(viridis)
library(RColorBrewer)
library(tidyr)
library(ggplot2)

setwd("E:/TIMELINE_SST/")

shp_path <- "GIS/Europe/Europe.gpkg"
plt_path = "OUT/Plots/Coast_distance/"
ocean_path <- "GIS/World_Seas_IHO_v3/"
europe_crs <- st_crs("+init=epsg:3035")

sites_path <- "GIS/Sites/"

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
