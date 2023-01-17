# Fetch cryo db data for the tabular version of the holocene arctic biodiversity
# map
# Jakob J. Assmann jakob.assmann@uzh.ch 17 January 2023

## Prepare envrionment ----
cat("Loading packages...")

# List of packages
required_packages <-
  c("terra")

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
locations <- read.csv("../visualisation/thalloo-static-site/map-data/ahbdb.txt")

test_rast <- rast("/vsicurl/https://os.zhdk.cloud.switch.ch/envicloud/chelsa/chelsa_V1/chelsa_trace/orog/CHELSA_TraCE21k_dem_-100_V1.0.tif")
