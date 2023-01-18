# Fetch cryo db data for the tabular version of the holocene arctic biodiversity
# map
# Jakob J. Assmann jakob.assmann@uzh.ch 17 January 2023

## Prepare envrionment ----
cat("Loading packages...")

# List of packages
library(terra)
library(sf)
library(pbapply)
library(tidyverse)

cat("done.\n")

# Work relative to cryo_db folder
# setwd("cryo_db") 

## Load site coordinates from ahbddb.txt
cat("Loading ahbdb coordinates...")
ahbdb <- read_delim("../visualisation/thalloo-static-site/map-data/ahbdb.txt",
                    delim = "\t",
                    show_col_types = FALSE)
# Parse and convert as vect
ahbdb <- st_as_sf(ahbdb, 
                  coords = c("LonDD", "LatDD"), 
                  crs = 4326,
                  remove = F) %>%
  select(LonDD, LatDD, geometry) %>%
  distinct()

cat("done.\n")

## Determine subset to extract
cat("Determining subest to extract...")

# Load already extracted coordinates
cryo_db <- read_csv("cryo_db.csv",
                    show_col_types = FALSE) %>%
  st_as_sf(coords = c("LonDD", "LatDD"), 
           crs = 4326,
           remove = F) 

# Determing subset of coordinates to extract
coords_to_extract <- ahbdb %>%
  filter(!(geometry %in% unique(cryo_db$geometry))) %>%
  # Add unique ID column
  mutate(ID = 1:nrow(.))

cat("done.\n")

# Status
cat("\nExtracting", nrow(coords_to_extract), "new records.",
    "\nTotal distinct locations in DB:", nrow(ahbdb),".\n\n")

## Load CHELSA TraCE21k data remotely

# Helper function to parse years into a more readable form
parse_year <- function(layer_names, var_name, new_var_name = NA){
  if(is.na(new_var_name)) new_var_name <- var_name
  new_layer_names <- layer_names %>% 
    gsub(paste0(".*", var_name, "_(.*)_.*"),
         paste0(new_var_name, "_\\1"), 
         .)
  return(new_layer_names)
}

# Mean temperature per century 
cat("Loading temperature data...")
CHELSA_temp <- read.csv("source_data/CHELSA_TraCE21k_bio01.txt") %>%
  lapply(function(x) paste0("/vsicurl/", x)) %>%
  lapply(rast) %>% 
  rast()
# Rename temperature layers
names(CHELSA_temp)<- names(CHELSA_temp) %>% parse_year("bio01", "temp")
cat("done.\n")  

# Mean precipitation per century
cat("Loading precipitation data...")
CHELSA_precip <- read.csv("source_data/CHELSA_TraCE21k_bio12.txt") %>%
  lapply(function(x) paste0("/vsicurl/", x)) %>%
  lapply(rast) %>% 
  rast()
# Rename precipitation layers
names(CHELSA_precip) <- names(CHELSA_precip) %>% parse_year("bio12", "precip")
cat("done.\n")  

# Digital elevation model
cat("Loading DEM data...")
CHELSA_dem <- read.csv("source_data/CHELSA_TraCE21k_dem.txt") %>%
  lapply(function(x) paste0("/vsicurl/", x)) %>%
  lapply(rast) %>% 
  rast()
# Rename dem layers
names(CHELSA_dem) <- names(CHELSA_dem) %>% parse_year("dem")
cat("done.\n")  

# Glacial extent
cat("Loading glacial extent polygons")
test <- st_read(paste0("/vsicurl/","https://chartercryodb.s3.eu-central-1.amazonaws.com/CHELSA_TraCE21_gle_geom/CHELSA_TraCE21k_gle_-101_V1.0.gpkg"))
CHELSA_gle_files <- list.files("source_data/CHELSA_TraCE21_gle_geom/",
                               full.names = T)
CHELSA_gle <- CHELSA_gle_files %>%
  lapply(read_sf) %>%
  setNames(gsub(".*/(CHELSA_TraCE21k_gle_.*).gpkg", "\\1", CHELSA_gle_files))
# Rename gle files
names(CHELSA_gle) <- names(CHELSA_gle) %>% parse_year("gle")
cat("done.\n")

## Load time-step look_up-table
time_id_look_up <- read.csv("source_data/time_id_look_up.csv")
# Helper function to look up k years BP
look_up_year <- function(time_id_vect){
  sapply(time_id_vect, 
         function(time_id){
           time_id_look_up$kBP[time_id == time_id_look_up$timeID]
         })
}

## Extract temperature, precipitation and elevation

# Temperature
cat("Extracting temperature values...")
croy_db_new_extract <- terra::extract(CHELSA_temp, coords_to_extract) %>% 
  pivot_longer(2:ncol(.), 
                              names_to = "year_kBP", 
                              values_to = "CHELSA_TraCE21k_temp") %>% 
  mutate(year_kBP = gsub("temp_", "", year_kBP)) %>%
  mutate(year_kBP = look_up_year(year_kBP))
cat("done.\n")  

# Precipitation
cat("Extracting precipitation values...")
croy_db_new_extract <- terra::extract(CHELSA_precip, coords_to_extract) %>% 
  pivot_longer(2:ncol(.), 
               names_to = "year_kBP", 
               values_to = "CHELSA_TraCE21k_precip") %>% 
  mutate(year_kBP = gsub("precip_", "", year_kBP)) %>%
  mutate(year_kBP = look_up_year(year_kBP)) %>%
  full_join(croy_db_new_extract, .)
cat("done.\n")  

# Elevation
cat("Extracting elevation values...")
croy_db_new_extract <- terra::extract(CHELSA_dem, coords_to_extract) %>% 
  pivot_longer(2:ncol(.), 
               names_to = "year_kBP", 
               values_to = "CHESA_TraCE21k_elevation") %>% 
  mutate(year_kBP = gsub("dem_", "", year_kBP)) %>%
  mutate(year_kBP = look_up_year(year_kBP)) %>%
  full_join(croy_db_new_extract, .)
cat("done.\n")  

# Calculate distance to land ice
cat("Calculating distance glaciers...")

# Helper function for distance calculation
croy_db_new_extract <- pblapply(seq_along(CHELSA_gle),
       function(index){
         # Convert coordinates to extract
         coords_to_extract_sf <- st_as_sf(coords_to_extract)
         # Get distances
         coords_to_extract_sf$CHESA_TraCE21k_dist_to_land_ice <- 
           st_distance(coords_to_extract_sf, CHELSA_gle[[index]]) %>%
           as.numeric()
         # Add year
         coords_to_extract_sf$year_kBP <- names(CHELSA_gle)[[index]]
         # Return as data frame
         coords_to_extract_sf %>%
           st_drop_geometry() %>%
           select(-LonDD, LatDD) %>%
           return()
       }) %>%
  reduce(full_join) %>%
  mutate(year_kBP = gsub("gle_", "", year_kBP)) %>%
  mutate(year_kBP = look_up_year(year_kBP)) %>%
  full_join(croy_db_new_extract, .)
cat("done.\n")

# # Calculate dates last glaciated and date last below sea-level
# croy_db_new_extract <- croy_db_new_extract %>%
#   group_by(id) %>%
  
## Add coordinates
croy_db_new_extract <- coords_to_extract %>%
  st_as_sf() %>%
  st_drop_geometry() %>%
  select(ID, LonDD, LatDD) %>%
  full_join(croy_db_new_extract)

## Merge with existing data
cryo_db <- bind_rows(croy_db_new_extract)

## Update database
write_csv(cryo_db, "cryo_db.csv")
