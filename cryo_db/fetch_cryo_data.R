# Fetch cryo db data for the tabular version of the holocene arctic biodiversity
# map
# Jakob J. Assmann jakob.assmann@uzh.ch 17 January 2023

## Prepare envrionment ----
cat("Loading packages...")

# List of packages
required_packages <-
  c("terra",
    "sf",
    "dplyr",
    "tidyr")

# Attempt load
packages_loaded <- sapply(required_packages, require, character.only = T)

# Check loading
if(all(packages_loaded)){
  cat("... all packages loaded.\n")
} else {
  # If not successful install and load packages that failed
  for(i in length(packages_loaded)){
    if(!packages_loaded[i]){
      cat("\nInstalling package", required_packages[i])
      install.packages(required_packages[i])
      library(required_packages[i])
      cat("Success.")
    }
  }
  cat("... all packages installed and loaded.\n")
}

## Load site coordinates from ahbddb.txt
ahbdb <- read.csv("../visualisation/thalloo-static-site/map-data/ahbdb.txt",
                      sep = "\t")
# Parse and convert as vect
ahbdb <- st_as_sf(ahbdb, 
                  coords = c("LonDD", "LatDD"), 
                  crs = 4326,
                  remove = F) %>%
  select(LonDD, LatDD, geometry) %>%
  distinct()

## Determine subset to extract
coords_to_extract <- ahbdb %>%
  mutate(ID = 1:nrow(.))
  
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

# Glatial extent
cat("Loading glacial extent data...")
CHELSA_gle <- read.csv("source_data/CHELSA_TraCE21k_gle.txt") %>%
  lapply(function(x) paste0("/vsicurl/", x)) %>%
  lapply(rast) %>% 
  rast()
# Rename gle layers
names(CHELSA_gle) <- names(CHELSA_gle) %>% parse_year("gle")
cat("done.\n")  

# # Glatial elevation
# cat("Loading glacier elevation data...")
# CHELSA_glz <- read.csv("source_data/CHELSA_TraCE21k_glz.txt") %>%
#   lapply(function(x) paste0("/vsicurl/", x)) %>%
#   lapply(rast) %>% 
#   rast()
# # Rename glz layers
# names(CHELSA_glz) <- names(CHELSA_gle) %>% parse_year("glt")
# cat("done.\n")  

## Load time-step look_up-table
time_id_look_up <- read.csv("source_data/time_id_look_up.csv")
# Helper function to look up k years BP
look_up_year <- function(time_id_vect){
  sapply(time_id_vect, 
         function(time_id){
           time_id_look_up$kBP[time_id == time_id_look_up$timeID]
         })
}

## Extract temperature, precipitation and eleveation

# Temperature
cat("Extracting temperature values...")
croy_db_new_extract <- terra::extract(CHELSA_temp, coords_to_extract[2,]) %>% 
  pivot_longer(2:ncol(.), 
                              names_to = "year_kBP", 
                              values_to = "CHELSA_TraCE21k_temp") %>% 
  mutate(year_kBP = gsub("temp_", "", year_kBP)) %>%
  mutate(year_kBP = look_up_year(year_kBP))
cat("done.\n")  

# Precipitation
cat("Extracting precipitation values...")
croy_db_new_extract <- terra::extract(CHELSA_precip, coords_to_extract[2,]) %>% 
  pivot_longer(2:ncol(.), 
               names_to = "year_kBP", 
               values_to = "CHELSA_TraCE21k_precip") %>% 
  mutate(year_kBP = gsub("precip_", "", year_kBP)) %>%
  mutate(year_kBP = look_up_year(year_kBP)) %>%
  full_join(croy_db_new_extract, .)
cat("done.\n")  

# Elevation
cat("Extracting elevation values...")
croy_db_new_extract <- terra::extract(CHELSA_dem, coords_to_extract[2,]) %>% 
  pivot_longer(2:ncol(.), 
               names_to = "year_kBP", 
               values_to = "CHESA_TraCE21k_elevation") %>% 
  mutate(year_kBP = gsub("dem_", "", year_kBP)) %>%
  mutate(year_kBP = look_up_year(year_kBP)) %>%
  full_join(croy_db_new_extract, .)
cat("done.\n")  

## Add coordinates
croy_db_new_extract <- coords_to_extract %>%
  st_as_sf() %>%
  st_drop_geometry() %>%
  select(ID, LonDD, LatDD) %>%
  full_join(croy_db_new_extract)

## Merge with existing data
cryo_db <- bind_rows(croy_db_new_extract)

## Update database
write.csv(cryo_db, "cryo_db.csv")
