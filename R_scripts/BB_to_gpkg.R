# write tile lists to geopackages

library(sf)

# path to tile lists
tile_path <- "E:/TIMELINE_SST/Tile_Lists/"
out_path <- "E:/TIMELINE_SST/GIS/bb_study_areas_3035.gpkg"

read_tile_lists <- function(tile_path, IHO_name){
  # get list of tile IDs covering IHO basin (from Drive)
  df <- read.csv(paste0(tile_path, IHO_name, '.csv'))
  return(df)
}

get_bb_list <- function(tile_df){
  # get bounding box from .csv files
  tile_bb <- st_bbox(c(xmin = min(tile_df$left), xmax = max(tile_df$right), 
                       ymax = max(tile_df$top), ymin = min(tile_df$bottom)), crs = 4326)
  tile_bb <- st_transform(st_as_sfc(tile_bb), crs = 3035)
  return(tile_bb)
}

poly1 <- read_tile_lists(tile_path, 'Adriatic Sea') %>% get_bb_list()
poly2 <- read_tile_lists(tile_path, 'Skagerrak') %>% get_bb_list()
poly3 <- read_tile_lists(tile_path, 'Aegean Sea') %>% get_bb_list()
poly4 <- read_tile_lists(tile_path, 'Balearic (Iberian Sea)') %>% get_bb_list()

st_write(c(poly1,poly2, poly3, poly4), out_path, layer = 'study_areas', driver = "GPKG")
