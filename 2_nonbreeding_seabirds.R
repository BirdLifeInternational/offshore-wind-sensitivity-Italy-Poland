#Use range maps geopackages to make seabird richness map
#Beth Clark May 2024 
rm(list=ls())

library(sf)
library(terra)

layer <- sf::st_layers("species_2024/species_2024.gpkg")$name[1]
layer

## read species list dataframe for non-breeding species
spp_list_df <- read.csv("nonbreeding_species.csv")

#Read in study raster
blank_ras <- terra::rast("seaward-ext-background-raster.tif")

#Non-breeding season ####

#native non-breeding is seasonal == 3
#(resident is seasonal == 1, native breeding is seasonal == 2)

#loop below took a while,  about 2 mins per species. 
for(i in 1:nrow(spp_list_df)){
  
  #select species
  query <- paste0(
    "SELECT * FROM \"", layer, "\" WHERE sci_name IN ('", spp_list_df$Scientific.name[i], "')"
  )
  
  #read in polygon for 1 species (takes a while)
  sp <- st_read('species_2024.gpkg', 
                layer = layer,
                query = query)
  
  #select the correct season
  sp <- subset(sp, seasonal == 3)
  
  sp <- st_transform(sp, crs = st_crs(blank_ras))
  
  #convert to raster
  sp_ras <- blank_ras
  sp_ras <- terra::mask(sp_ras, vect(sp))
  plot(sp_ras, main = spp_list_df$Scientific.name[i])
  plot(eez[1], add = T, color = "transparent")
  
  #save raster
  writeRaster(sp_ras, paste0(folder, 
                             spp_list_df$Scientific.name[i],
                             ".tif"), overwrite = T)
  
  #Bradbury method to weight species by sensitivity 
  #(max score of collision and displacement)
  species_sensitivity_raster <- log(sp_ras+1)*
    max(spp_list_df$collision_score[spp_list_df$species_label == species_label],
        spp_list_df$displacement_score[i][spp_list_df$species_label == species_label],
        na.rm = T) 

  #combine collision
  if(i == 1){
    all_species <- species_sensitivity_raster
  } else {
    all_species <- sum(c(all_species,species_sensitivity_raster), na.rm = T)
  }
  
  print(i)
}


plot(all_species)

writeRaster(all_species, "Birdlife_range_MaxComb_Disp_Coll.tif", overwrite = T)
