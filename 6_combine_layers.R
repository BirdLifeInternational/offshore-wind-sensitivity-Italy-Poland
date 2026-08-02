#Beth Clark 2023
#Convert seaward extension rasters to sensitivity maps and combine using Bradbury method
#The combine with additional layers
rm(list=ls()) 

library(sf)
library(terra)
library(tidyverse)
library(RColorBrewer)
library(BAMMtools)

land <- read_sf("land")
eez <- read_sf("EEZ_land_union_v3_202003")

tif_folder <- "seaward_extension_outputs\\"
rasters <- list.files(tif_folder, pattern = ".tif");rasters

species_raster <- rast(paste0(tif_folder,rasters[1]))

#Transform shapefiles for plotting
it_eez <- st_transform(it_eez, crs = st_crs(species_raster))
land <- st_transform(land, crs = st_crs(species_raster))
land_it <- st_transform(land_it, crs = st_crs(species_raster))

#Create colour ramps
greens <- colorRampPalette(brewer.pal(n = 9, name = "Greens"))(100)[15:100]

#Read in sensitivity scores
#csv containing columns named "Scientific_name", "Collision_SI", "Displacement_SI"
sensitivity <- read.csv("sensitivity.csv")
head(sensitivity)

i<-1

#start from 2 as 1 is all spp preening buffer
for(i in 2:length(rasters)){
  
  species_raster <- rast(paste0(tif_folder,rasters[i]))
  
  #Extract species name from file name
  species <- substring(sources(species_raster),134,nchar(sources(species_raster))-4);species
  
  #Trim raster to EEZ
  species_raster <- terra::mask(species_raster, vect(it_eez))
  
  #Normalise each map to max 1
  species_raster <- species_raster/max(values(species_raster),na.rm = T)
  
  #Bradbury method to weight species by sensitivity 
  #(max score of collision and displacement)
  species_sensitivity_raster <- log(species_raster+1)*
    max(sensitivity$Collision_SI[sensitivity$species_label == species_label],
        sensitivity$Displacement_SI[sensitivity$species_label == species_label],
        na.rm = T)    
  
  #Save rasters 
  if(exists("species_sensitivity_raster_stack")){
    species_sensitivity_raster_stack <- c(species_sensitivity_raster_stack,
                                          species_sensitivity_raster)
  } else {
    species_sensitivity_raster_stack <- species_sensitivity_raster
  }
  
  print(i)
  print(species)
}

#combine into signal map using Bradbury formula
breeding_map <- sum(species_sensitivity_raster_stack,
                    na.rm = T)
summary(values(breeding_map))

#Normalise to max 1
seabirds_breeding <- breeding_map/max(values(breeding_map),na.rm = T)

avistep <- c("white",colorRampPalette(c("#55FF00", "#FFFF00", "#FFAA00", "#FF0000"))(1000))
avistep4 <- colorRampPalette(c("#55FF00", "#FFFF00", "#FFAA00", "#FF0000"))(4)

#quick plot to check
plot(seabirds_breeding, col = avistep)
plot(it_eez[1], add = T, col = "transparent", border = "blue")
plot(land, add = T, col = "grey90", border = "grey85")

#combine with non-breeding seabirds layer from BirdLife range maps ####
bl_ranges <- rast("Birdlife_range_MaxComb_Disp_Coll.tif")
plot(bl_ranges)
length(unique(values(bl_ranges)))

#crop to match
seabirds_breeding <- crop(seabirds_breeding,bl_ranges)
length(unique(values(seabirds_breeding)))

#plot
par(mfrow = c(1,2))
plot(seabirds_breeding, col = avistep)
plot(land, add = T, col = "grey90", border = "grey85")

plot(bl_ranges, col = avistep)
plot(land, add = T, col = "grey90", border = "grey85")

#combine with max of all layers
br_nonbr_seabirds <- max(seabirds_breeding,bl_ranges)

plot(br_nonbr_seabirds, col = avistep)
plot(land, add = T, col = "grey90", border = "grey85")

#add in 5km preening/rafting buffer for all seabird nesting sites
seabird_colonies <- rast("seaward-ext-multi-col-dist-all-species-5km.tif")
seabird_cols <- crop(seabird_colonies, seabirds_breeding)
seabird_cols[is.na(seabird_cols)] <- 0
plot(seabird_cols)

#combine with max of all layers
all_seabirds <- max(seabirds_breeding,bl_ranges,seabird_cols)

plot(all_seabirds, col = avistep)
plot(land, add = T, col = "grey90", border = "grey85")

par(mfrow = c(1,2))
plot(all_seabirds, col = avistep)
plot(land, add = T, col = "grey90", border = "grey85")

#Next step is to add other layers ####

#read in raster for migratory large landbird (from ISPRA 2021)
mig_large_landbirds <- rast("mig_large_landbirds.tif")

#20km+10km buffer around small islands (from ISPRA 2021)
small_islands <- rast("Small_Islands.tif")
plot(small_islands)

#read in IBA polygons
ibas_new <- st_read("IBAs")
ibas_new
plot(ibas_new[1])
ibas_new <- st_transform(ibas_new, crs = st_crs(species_raster))

#rasterise to match species data
ibas_mask <- species_raster
values(ibas_mask) <- 1 
vect(ibas_new[1])
ibas_mask <- terra::mask(ibas_mask, vect(ibas_new[1]))
plot(ibas_mask, col = "darkred")

#combine with seabird layers
all_layers <- max(all_seabirds,mig_large_landbirds,small_islands,ibas_mask)

#apply Jenk's natural breaks ####
breaks <- getJenksBreaks(values(all_layers), k = 5);breaks

all_offshore_natural_breaks <- all_layers
values(all_offshore_natural_breaks) <- ifelse(values(all_layers) <= breaks[2], 1, values(all_offshore_natural_breaks)) 
values(all_offshore_natural_breaks) <- ifelse(values(all_layers) > breaks[2] & values(all_layers) <= breaks[3], 2, values(all_offshore_natural_breaks)) 
values(all_offshore_natural_breaks) <- ifelse(values(all_layers) > breaks[3] & values(all_layers) <= breaks[4], 3, values(all_offshore_natural_breaks)) 
values(all_offshore_natural_breaks) <- ifelse(values(all_layers) > breaks[4], 4, values(all_offshore_natural_breaks)) 

par(mfrow = c(1,1))
plot(all_offshore_natural_breaks, col = avistep4)
plot(it_eez[1], add = T, col = "transparent", border = "blue")
plot(land, add = T, col = "grey90", border = "grey85")
