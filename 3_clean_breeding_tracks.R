#Bethany Clark, BirdLife International
rm(list=ls())

#load libraries ####
library(tidyverse)
## sf package for spatial data analyses (i.e. vector files such as points, lines, polygons)
library(sf)
## Tidyverse for data manipulation & plots
library(tidyverse)
## rnaturalearth package for geographic basemaps in R
library(rnaturalearth);library(rnaturalearthdata)
## lubridate for date time
library(lubridate)
## speed filter
library(trip)
## linear interpolation
library(adehabitatLT)
library(track2KBA)

#read in high resolution large scale of relevant countries: the 
#study country and surrounding countries used by tracked birds
## Land polygon from https://gadm.org/
basemap <- read_sf("gadm_410.gpkg")
unique(basemap$NAME_0)
countries <- c("Country1","Country2")
land <- subset(basemap, NAME_0 %in% countries)

#Read in Exclusive Economic Zones file
eez_file <- read_sf("EEZ_land_union_v3_202003")
eez_it <- subset(eez_file, SOVEREIGN1 == "Country")
it_extent <- st_bbox(eez_it)

#make sure plot includes entire IT EEZ and all colonies
xmin <- it_extent$xmin[[1]]
xmax <- it_extent$xmax[[1]]
ymin <- it_extent$ymin[[1]]
ymax <- it_extent$ymax[[1]]
border <- 0.01

plot(eez_it[1])

dir_in <- "All_tracking_STDB_format/"
dir.create("./New_tracking_data_analysis/1_maps/")
dir.create("./New_tracking_data_analysis/1_clean_tracks_per_dataset/")
dir.create("./New_tracking_data_analysis/2_interpolated_tracks_per_dataset/")

#read files ####
files <- list.files(dir_in);files
meta_data <- read.csv("tracking_meta_data.csv")
tail(meta_data)

species <- unique(meta_data$scientific_name);species

meta_data$species_colony <- paste(meta_data$common_name, meta_data$colony_name, sep="_")

meta_data <- meta_data %>% 
  arrange(species_colony) 

#add population info
colonies_meta <- as.data.frame(files)
colonies_meta$n_birds <- NA
colonies_meta$n_trips <- NA
colonies_meta$col_dist_max_km <- NA
colonies_meta$interp_interval <- NA
colonies_meta$interp_interval_max <- NA

#choose colony ####
length(files)

for(j in c(1:length(files))){ #change for elonora's
  
  #read in files in the format downloaded from seabirdtracking.org
  df_stdb_output <- read.csv(paste0(dir_in,files[j]))
  head(df_stdb_output)
  
  table(df_stdb_output$dataset_id)
  
  #Subset only breeding adults
  df_stdb_output <- subset(df_stdb_output, age == "adult")
  df_stdb_output <- subset(df_stdb_output, breed_stage != "non-breeding")
  df_stdb_output <- subset(df_stdb_output, breed_stage != "fail (breeding season)")
  df_stdb_output <- subset(df_stdb_output, breed_stage != "pre-egg")
  df_stdb_output <- subset(df_stdb_output, breed_stage != "unknown") 
  
  table(df_stdb_output$dataset_id)
  
  head(df_stdb_output)
  
  #fix datetime if in old format
  df_stdb_output$date_gmt <- ifelse(grepl("/",df_stdb_output$date_gmt),
                                    paste(substr(df_stdb_output$date_gmt,7,10),
                                          substr(df_stdb_output$date_gmt,4,5),
                                          substr(df_stdb_output$date_gmt,1,2),
                                          sep="-"),
                                    df_stdb_output$date_gmt)
  
  df_stdb_output$dttm <- with(df_stdb_output, ymd(date_gmt) + hms(time_gmt))
  
  ## first check how many duplicate entries you may have. If there are many, it
  ## is worth exploring your data further to understand why.
  n_duplicates <- df_stdb_output %>% 
    group_by(bird_id) %>% 
    arrange(dttm) %>% 
    dplyr::filter(duplicated(dttm) == T)
  
  ## review how many duplicate entries you may have. Print the message:
  print(paste("you have ",nrow(n_duplicates), " duplicate records in a dataset of ",
              nrow(df_stdb_output), " records.", sep =""))
  rm(n_duplicates)
  
  ## remove duplicates entries if no further exploration is deemed necessary
  df_stdb_output_no_dups <- df_stdb_output %>% 
    ## first group data by individual animals and unique track_ids
    group_by(bird_id, track_id) %>% 
    ## then arrange by timestamp
    arrange(dttm) %>% 
    ## then if a timestamp is duplicated (TRUE), then don't select this data entry.
    ## only select entries where timestamps are not duplicated (i.e. FALSE)
    dplyr::filter(duplicated(dttm) == F)
  
  table(df_stdb_output$bird_id)
  table(df_stdb_output_no_dups$bird_id)

  df_stdb_output <- df_stdb_output_no_dups
  
  #save out plot in a folder for checking
  sp_dat <- df_stdb_output %>%
    st_as_sf(coords = c("longitude","latitude"), crs = 4326) 
  png(paste0("./New_tracking_data_analysis/1_maps/",j,".png"),
      width = 900, height = 850,)
  ggplot(eez_it) +
    geom_sf(alpha = 0.1, fill = "blue")+
    geom_sf(data = land, fill = "darkgrey")+
    geom_sf(data = sp_dat, 
            alpha = 0.6, size = 1.2, shape = 20)+
    ggtitle(files[j])
  dev.off()
  
  ## Define maximum speed in km/h (kilometers per hour)
  speed.filter.threshold <- 100 ## Important number to change depending on whether you have flying or non-flying seabirds
  
  ## create blank data frame to capture filtered tracks and summary data
  tracks_speed <- data.frame()
  tracks_speed_summary <- data.frame()

  ## 
  for(i in 1:length(unique(df_stdb_output$bird_id))){
    temp <- df_stdb_output %>% dplyr::filter(bird_id == unique(df_stdb_output$bird_id)[i])
    
    #only use if 5 or more locations
    if(nrow(temp) > 4){
      
      ## remove any erroneous locations due to speed use the McConnel Speed Filter 
      ##from the trip package
      trip_obj <- temp %>% 
        group_by(bird_id) %>% 
        dplyr::select(x = latitude, 
                      y = longitude, 
                      dttm, 
                      everything()) %>% 
        trip()
      
      ## McConnel Speedilter -----
      ## apply speedfilter and create data frame
      trip_obj$Filter <- speedfilter(trip_obj, max.speed = speed.filter.threshold)  # speed in km/h
      trip_obj <- data.frame(trip_obj)
      #head(trip_obj,2)
      #dim(trip_obj)
      
      ## Keep only filtered coordinates - after checking dimensions of other outputs again
      trip_obj <- subset(trip_obj,trip_obj$Filter==TRUE)
      
      ## bind back onto dataframe
      tracks_speed <- rbind(tracks_speed, trip_obj)
      
      ## Populate summary data
      temp_summary <- data.frame(bird_id = temp$bird_id[1],
                                 n_points_PreFilter = nrow(temp),
                                 n_points_PostFilter = nrow(trip_obj),
                                 points_removed = ifelse(nrow(temp)-nrow(trip_obj) > 0, 
                                                         nrow(temp)-nrow(trip_obj), "-"))
      
      ## bind on summary information
      tracks_speed_summary <- rbind(tracks_speed_summary, temp_summary)
      
      ## remove temporary items before next loop iteration
      rm(temp,trip_obj, temp_summary)
    } else {
      
      ## Populate summary data
      temp_summary <- data.frame(bird_id = temp$bird_id[1],
                                 n_points_PreFilter = nrow(temp),
                                 n_points_PostFilter = 0,
                                 points_removed = nrow(temp))
      
      ## bind on summary information
      tracks_speed_summary <- rbind(tracks_speed_summary, temp_summary)
    }
    
    ## Print loop progress
    print(paste("Track ", i, " of ", length(unique(df_stdb_output$bird_id)), " processed"))
    
  }
  #warnings checked - they are fine
  head(tracks_speed)
  tracks_speed_summary
  
  ## update column names in speed filtered tracks
  tracks_speed <- tracks_speed %>% 
    mutate(latitude = x,
           longitude = y)
  
  ## Format the key data fields to the standard used in track2KBA
  dataGroup <- formatFields(
    ## your input data.frame or tibble
    dataGroup = tracks_speed, 
    ## ID of the animal you tracked
    fieldID   = "bird_id", 
    ## date time in GMT
    fieldDateTime = "dttm", 
    ## longitude of device
    fieldLon  = "longitude", 
    ## latitude of device
    fieldLat  = "latitude"
  )
  dataGroup$dttm <- dataGroup$DateTime
  
  ## Check output. Output is a data.frame
  head(dataGroup,2)
  
  #trip split ####
  colony <- data.frame(Longitude = dataGroup$lon_colony[1], 
                       Latitude  = dataGroup$lat_colony[1])
  
  ## Check output. Output is a data.frame
  head(dataGroup,2)
  str(dataGroup)
  
  ## First define your key parameters outside of the function. Useful for using them later again if needed.
  inner.buff.distance = 3 # km - defines distance an animal must travel to count as trip started
  return.buff.distance = 10 # km - defines distance an animal must be from the colony to have returned and thus completed a trip
  duration.time = 1 # hours - defines time an animal must have traveled away from the colony to count as a trip. helps remove glitches in data or very short trips that were likely not foraging trips.
  
  ## Input is a 'data.frame' of tracking data and the central-place location(s). 
  ## Output is a 'SpatialPointsDataFrame'.
  trips <- tripSplit(
    dataGroup  = dataGroup,
    colony     = colony, # define source location.
    innerBuff  = inner.buff.distance,      
    returnBuff = return.buff.distance,     
    duration   = duration.time,     
    nests = F,           # specify nests = T if using unique colony locations per animal,
    gapLimit = NULL, # The period of time between points (in days) to be considered too large to be a contiguous tracking event
    rmNonTrip  = T    # If true, points not associated with a trip will be removed / if false, points not associated with a trip will be kept
  )
  
  #toollkit example using rmNonTrip  = F, but this doesn't work if you don't want to remove partial trips
  
  table(trips$tripID)
  
  ## Review data after tripSplit()
  head(trips,2)
  
  table(trips$Returns)
  
  ## Split the points
  trips2 <- trips[trips$ColDist > inner.buff.distance*1000, ]
  plot(trips$X,trips$Y)
  plot(trips2$X,trips2$Y)
  
  mapTrips(trips = trips2, colony = colony)
  
  trips_df <- as.data.frame(trips2)
  trips_df$bird_id <- trips_df$ID
  
  head(trips_df)
  
  ## summary of trips associated with Return or not
  totalTripsAll <- data.frame(trips) %>% group_by(tripID, Returns) %>% 
    summarise(count = n()) %>% 
    data.frame(.)
  ## view summary result
  table(totalTripsAll$Returns)
  
  trips <- trips_df
  
  ## create data frame and find IDs of trips with >5 locations; as required for track2KBA analysis
  trips_to_keep <- data.frame(trips) %>% 
    group_by(tripID) %>% 
    summarise(triplocs = n()) %>% 
    dplyr::filter(triplocs > 5)
  
  ## select the relevant tripIDs only and create new object
  trips_df <- data.frame(trips) %>% 
    dplyr::filter(tripID %in% trips_to_keep$tripID)
  
  
  ## tripSummary() ----
  sumTrips <- tripSummary(trips = trips_df, colony = colony, nests = T)
  
  ## Compare how many trips you removed
  ## Before
  length(unique(trips$tripID))
  
  ## After
  length(unique(trips_df$tripID))   
  
  #save clean tracks
  
  write.csv(trips_df,paste0(
                "New_tracking_data_analysis/1_clean_tracks_per_dataset/",
                gsub(".csv","_",files[j]), "cleaned.csv"), row.names = F)
  
  colonies_meta$n_birds[j] <- length(unique(trips_df$tripID)) 
  colonies_meta$n_trips[j] <- length(unique(trips_df$ID)) 
  colonies_meta$col_dist_max_km[j] <- max(trips_df$ColDist)/1000
  
  #time diffs ####
  timeDiff <- trips_df %>% 
    data.frame() %>% 
    group_by(tripID) %>% 
    arrange(DateTime) %>% 
    mutate(delta_secs = as.numeric(difftime(DateTime, 
                                            lag(DateTime, default = first(DateTime)), 
                                            units = "secs"))) %>% 
    slice(2:n()) 
  
  head(data.frame(timeDiff))
  
  ## Summarise results by tripID
  SummaryTimeDiff <- timeDiff %>% 
    group_by(ID) %>% 
    summarise(mean_timegap_secs = mean(delta_secs),
              median_timegap_secs = median(delta_secs),
              min_timegap_secs = min(delta_secs),
              max_timegap_secs = max(delta_secs)) %>%
    ## time in days
    mutate(max_timegap_days =  max_timegap_secs / 86400) %>% 
    mutate(max_timegap_days = round(max_timegap_days,2)) %>% 
    data.frame()
  
  ## View results
  SummaryTimeDiff
  
  ## Average sampling interval of all data, median of median time gaps in minutes
  
  #interp interval ####
  interp.interval <- median(SummaryTimeDiff$median_timegap_secs);interp.interval
  colonies_meta$interp_interval[j] <- interp.interval
  colonies_meta$interp_interval_max[j] <- max(SummaryTimeDiff$median_timegap_secs)
  
  ## 
  trips_meta <- trips_df %>% 
    dplyr::select(scientific_name,
                  common_name,
                  site_name,
                  colony_name,
                  lat_colony,
                  lon_colony,
                  bird_id = ID,
                  trip_id = tripID,
                  age,
                  sex,
                  breed_stage,
                  breed_status) %>% 
    group_by(trip_id) %>% 
    slice(1)
  
  trips_interp_df <- data.frame()
  
  head(trips_df,2)
  
  ## Linear interpolation -----
  for(i in 1:length(unique(trips_df$tripID))){
    
    temp <- trips_df %>% dplyr::filter(tripID == unique(trips_df$tripID)[i])
    
    ## Apply linear interpolation step to speed filtered only data
    
    ## create ltraj object
    trip_lt <- as.ltraj(xy = bind_cols(x = temp$Longitude, 
                                       y = temp$Latitude),
                        date = temp$DateTime ,
                        id = temp$tripID)
    
    ## Linearly interpolate/re-sample tracks according to chose interval (specified in seconds)
    trip_interp <- redisltraj(trip_lt, 
                              interp.interval, 
                              type="time")
    head(trip_interp)
    
    ## convert back into format for track2KBA - dataframe for now
    trip_interp <- ld(trip_interp) %>% 
      dplyr::mutate(Longitude = x,
                    Latitude = y)
    
    ## bind back onto dataframe
    trips_interp_df <- rbind(trips_interp_df, trip_interp)
    
    ## remove temporary items before next loop iteration
    rm(temp,trip_lt)
    
    ## Print loop progress
    print(paste("Trip ", i, " of ", length(unique(trips_df$tripID)), " processed"))
    
  }
  
  ## First create columns with same names
  trips_interp_df <- trips_interp_df %>% 
    rename(trip_id = id)
  
  trips_meta$trip_id <- as.character(trips_meta$trip_id)
  
  ## Now bind the metadata back
  trips_interp_df <- left_join(trips_interp_df, 
                               trips_meta,
                               by = "trip_id")
  
  ## Review data
  head(trips_interp_df,2)
  
  ## Compare and ensure same number of trips across data
  length(unique(trips_interp_df$trip_id))
  
  #Save csv
  
  write.csv(trips_df, paste0(
    "New_tracking_data_analysis/2_interpolated_tracks_per_dataset/",
    gsub(".csv","_",files[j]), "interp.csv"), row.names = F)
  
  ## Format the data BACK INTO key data fields to the standard used in track2KBA
  dataGroup_interp <- formatFields(
    ## your input data.frame or tibble
    dataGroup = trips_interp_df, 
    ## ID of the animal you tracked
    fieldID   = "bird_id", 
    fieldDateTime = "date",
    ## longitude of device
    fieldLon  = "Longitude", 
    ## latitude of device
    fieldLat  = "Latitude"
  )
  
  ## Check output. Output is a data.frame
  head(dataGroup_interp,2)
  
  write.csv(colonies_meta,"./New_tracking_data_analysis/Colonies_metadata.csv",
            row.names = F)
  
  print(j)
}    
