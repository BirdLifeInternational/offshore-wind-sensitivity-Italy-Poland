#Bethany Clark, BirdLife International
rm(list=ls())

#load packages ####
library(tidyverse)
## sf package for spatial data analyses (i.e. vector files such as points, lines, polygons)
library(sf)
## Tidyverse for data manipulation & plots
library(tidyverse)
## leaflet package for interactive maps in R
library(leaflet)
## lubridate for date time
library(lubridate)
## speed filter
library(trip)
## linear interpolation
library(adehabitatLT)
##kde
library(adehabitatHR)
#rasters
library(terra)
library(track2KBA)
#colours for plots
library(RColorBrewer)
greens <- c("white",colorRampPalette(brewer.pal(n = 9, name = "Greens"))(100)[15:100])

#read in project grid
base_grid <- terra::rast("seaward-ext-background-raster.tif")
grid_cell_size <- 5000

#make a grid based on base_grid in the correct format
c <- as.numeric(ext(base_grid)[1]+2500)  ## to check my min x
d <- as.numeric(ext(base_grid)[2])  ## to check my max xx

e <- as.numeric(ext(base_grid)[3]+2500)   ## to check my min y
f <- as.numeric(ext(base_grid)[4])   ## to check my max y

a <- seq(c, d, by = grid_cell_size)
b <- seq(e, f, by = grid_cell_size)
null_grid <- expand.grid(x=a,y=b)
sp::coordinates(null_grid) <- ~x+y
sp::gridded(null_grid) <- TRUE
class(null_grid)

land <- read_sf("land")
eez <- read_sf("EEZ_land_union_v3_202003")

eez_file <- read_sf("input_data/EEZ_land_union_v3_202003")
eez_it <- subset(eez_file, SOVEREIGN1 == "Country")
it_extent <- st_bbox(eez_it)
#make sure plot includes entire EEZ and all colonies
xmin <- it_extent$xmin[[1]]
xmax <- it_extent$xmax[[1]]
ymin <- it_extent$ymin[[1]]
ymax <- it_extent$ymax[[1]]
border <- 0.01

plot(eez_it[1])

dir_in <- "./2_interpolated_tracks_per_dataset/"
dir.create("./3_kde_rasters/")
dir.create("./3_kde_maps/")

#read files ####
files <- list.files("2_interpolated_tracks_per_dataset/");files
meta_data <- read.csv("tracking_meta_data.csv")

#meta data files contains colony population sizes
head(meta_data)
tail(meta_data)

species <- unique(meta_data$scientific_name);species

meta_data$species_colony <- paste(meta_data$common_name, meta_data$colony_name, sep="_")

meta_data <- meta_data %>% 
  arrange(species_colony) 

colonies <- unique(meta_data$species_colony)

#create a dataframe to store results
colonies_meta <- as.data.frame(colonies)
colonies_meta$col_size <- NA
colonies_meta$n_files <- NA
colonies_meta$files <- NA
colonies_meta$n_birds <- NA
colonies_meta$n_trips <- NA
colonies_meta$col_dist_max_km <- NA
colonies_meta$interp_interval <- NA
colonies_meta$interp_interval_max <- NA
colonies_meta$h_ref <- NA
colonies_meta$h <- NA

#loop through all colonies####

for(j in c(1:length(colonies))){
  sp_col_df <- subset(meta_data, species_colony2 == colonies[j]);sp_col_df
  
  print(sp_col_df)
  
  colonies_meta$n_files[j] <- nrow(sp_col_df)
  colonies_meta$files[j] <- paste(sp_col_df$files, collapse = "_")
  
  colonies_meta$col_size[j] <- sp_col_df$pop_size[1]
  
  #read all the tracking data files for the colony and rbind them
  trips_interp_df<-do.call("rbind",lapply(as.character(paste0(dir_in, 
                                                              sp_col_df$files_interp)),
                                          read.csv,stringsAsFactors = F)) 
  
  plot(trips_interp_df$Longitude,trips_interp_df$Latitude)
  
  trips_interp_df$dttm <- ymd_hms(trips_interp_df$dttm)
  trips_interp_df <- trips_interp_df %>%
    drop_na(dttm)
  
  #reproject to match the base grid ####
  sp_interp <- trips_interp_df %>%
    st_as_sf(coords = c("Longitude","Latitude"),
             crs = 4326) 
  dat_proj <- st_transform(sp_interp, crs = st_crs(base_grid))
  
  dat_points <- as.data.frame(st_coordinates(dat_proj))
  
  xy <- SpatialPoints(dat_points)
  
  crs(xy) <- crs(base_grid)
  
  #plot(xy)
  
  #Create UD with KDE ####
  kud <- kernelUD(xy, grid = null_grid)
  kud@h
  colonies_meta$h_ref[j] <- kud@h$h
  # h = href is the default - ad hoc method for determining h
  
  write.csv(colonies_meta,"./Colonies_metadata.csv",
            row.names = F)
  
  print(j)
}    

#numbers below need to be updated if running again after adding eleonora's
#calculate the mean href for the species and redo.
species1_h <- round(mean(colonies_meta$h_ref[1:9]),-2);species1_h
species2_h <- round(mean(colonies_meta$h_ref[10]),-2);species2_h

j = 11
#calc final kde rasters ####
for(j in c(10:12)){
  sp_col_df <- subset(meta_data, species_colony2 == colonies[j]);sp_col_df
  
  colonies_meta$col_size[j] <- sp_col_df$pop_size[1]
  
  print(sp_col_df)
  
  trips_interp_df<-do.call("rbind",lapply(as.character(paste0(dir_in, 
                                                              sp_col_df$files_interp)),
                                          read.csv,stringsAsFactors = F)) 
  
  plot(trips_interp_df$Longitude,trips_interp_df$Latitude)
  
  trips_interp_df$dttm <- ymd_hms(trips_interp_df$dttm)
  trips_interp_df <- trips_interp_df %>%
    drop_na(dttm)
  
  #reproject to match the base grid ####
  sp_interp <- trips_interp_df %>%
    st_as_sf(coords = c("Longitude","Latitude"),
             crs = 4326) 
  dat_proj <- st_transform(sp_interp, crs = st_crs(base_grid))
  
  dat_points <- as.data.frame(st_coordinates(dat_proj))
  
  xy <- SpatialPoints(dat_points)
  
  crs(xy) <- crs(base_grid)
  
  plot(xy)
  
  #Create UD with KDE ####
  
  #choose the species specific h (smoothing parameter)
  if(sp_col_df$scientific_name[1] == "Species1"){
    h_value <- species1_h
  } else if(sp_col_df$scientific_name[1] == "Species2"){
    h_value <- species2_h
  }
  
  colonies_meta$h[j] <- h_value
  
  kud <- kernelUD(xy, grid = null_grid, h = h_value)

  kde_spixdf <- kud[1]
  stk <- raster::stack(kde_spixdf)
  
  rast <- stk[[1]]
  rast_pop <- rast/max(values(rast))*colonies_meta$col_size[j]
  plot(rast_pop)
  
  rast_name <- paste0("./3_kde_rasters/",j,"_",
                      gsub("/","_",colonies[j]),
                      ".tif")
  
  raster::writeRaster(rast_pop, filename = rast_name, 
                      format = "GTiff", overwrite = TRUE)
  
  write.csv(colonies_meta,"./Colonies_metadata.csv",
            row.names = F)
  
  plot(rast_pop, col = greens, main = species)
  plot(eez[1], add = T, col = "transparent", border = "blue")
  plot(land, add = T, col = "grey", border = "grey")
  plot(land_it[1], add = T, col = "darkgrey", border = "black")  
  
  print(j)
}    
