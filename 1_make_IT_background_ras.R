#Bethany Clark, BirdLife International
rm(list=ls())

library(terra)
library(sf)

#Ready in study area polygon in suitable equal areas project (unit = m)
grid_shp <- read_sf("study_area_polygon") 

bounds <- st_bbox(grid_shp)

## Create Raster
ras <- terra::rast(grid_shp)

rm(grid_shp)

ras <- terra::rast(xmin = bounds[[1]], 
                   ymin = bounds[[2]], 
                   xmax = bounds[[3]], 
                   ymax = bounds[[4]], 
                   resolution = 5000,
                   crs = st_crs(ras)$wkt)
values(ras) <- 0

## Land polygon from https://gadm.org/
basemap <- read_sf("gadm_410.gpkg")

unique(basemap$NAME_0)

countries <- c("Country1","Country2")

basemap2 <- subset(basemap, NAME_0 %in% countries)

#If needed, convert land map to basemap projection
if(st_crs(basemap2) != st_crs(ras)){
  land <- st_transform(basemap2, crs = st_crs(ras))
} else {
  land <- basemap2
}

land2 <- st_union(land)

st_write(land2,"land.shp")

basemap_vector <- vect(land)
mask <- terra::rasterize(land, ras)
mask[is.na(mask)] <- 0
mask[mask == 1] <- NA
mask[mask == 0] <- 1
plot(mask)

terra::writeRaster(mask,"seaward-ext-background-raster.tif", 
                    overwrite = TRUE)
