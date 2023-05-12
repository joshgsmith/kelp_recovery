#Joshua G. Smith 
#jossmith@mbayaq.org


# Read data
################################################################################

# Clear workspace
rm(list = ls())

librarian::shelf(rerddap, lubridate, tidyverse, beepr, data.table, purrr, sf)


# Directories
basedir <- "/Volumes/seaotterdb$/kelp_recovery"
figdir <- here::here("analyses","4patch_drivers","Figures")

# Source
# https://coastwatch.pfeg.noaa.gov/erddap/

# Multi-scale Ultra-high Resolution (MUR) SST Analysis fv04.1, Global, 0.01°, 2002-present, Daily

# Build data
################################################################################

# Look up datasets 
datasets_grids <- ed_datasets(which="griddap", url="https://coastwatch.pfeg.noaa.gov/erddap/")
datasets_tables <- ed_datasets(which="tabledap", url="https://coastwatch.pfeg.noaa.gov/erddap/")


# Get data
################################################################################

#set export location
export_path <- file.path("SET_EXPORT_PATH")

#Set start and end dates
start_date <- as.Date("2003-01-01")
end_date <- as.Date("2020-12-31")


# Loop over years
for (year in unique(year(start_date):year(end_date))) {
  # Calculate start and end dates for current year
  chunk_start <- as.Date(paste0(year, "-01-01"))
  chunk_end <- as.Date(paste0(year, "-12-31"))
  
  # Download daily sst data for current year
  data_info <- info("jplMURSST41")
  data_chunk <- griddap(dataset="jplMURSST41", 
                        time = c(chunk_start, chunk_end),
                        longitude = c(-125, -117), latitude = c(36.500888, 37.122990), #Ano Nuevo to Lobos and offshore 5 Km
                        fields="analysed_sst")
  
  # Combine data frames in data_list
  data_list <- list(data_chunk)
  data_combined <- do.call(rbind, lapply(data_list, function(x) x$data))
  
  # Remove rows with missing values
  data_combined <- data_combined[complete.cases(data_combined), ]
  
  #Make into df and clean
  sst_daily <- data_combined %>%
    mutate(date=gsub("T00:00:00Z", "", time),
           year = lubridate::year(date),
           month = lubridate::month(date),
           day = lubridate::day(date))
  
  # Create a filename for current year
  filename <- paste0(year, "_sst_daily.Rds")
  
  # Save data for current year to an Rds file
  saveRDS(sst_daily, file.path(export_path,filename))
  
  # Remove current run from local environment to save memory 
  rm(data_info)
  rm(data_list)
  rm(sst_daily)
  rm(data_combined)
  rm(data_chunk)
}


# Read .Rds files inside directory and merge
################################################################################

setwd("/Users/jossmith/Desktop/test_data")

# Create a list of all the .Rds files in the directory
file_list <- list.files(pattern = "*.Rds")

# Read all the .Rds files into a list using purrr::map()
file_list <- map(file_list, readRDS)

# Merge all the .Rds files using data.table::rbindlist()
merged_data <- rbindlist(file_list)

# Save the merged data as a single .Rds file
saveRDS(merged_data, "2003_2022_daily_SST.Rds")


# Test load to local environment
################################################################################

sst_dat <- readRDS("2003_2022_daily_SST.Rds")

names(sst_dat)



# calculate SST at central coast pisco sites
################################################################################
#intersect monitoring sites with SST

#load sites --- available from https://opc.dataone.org/view/MLPA_kelpforest.metadata.2
site_table <- read.csv(file.path(basedir, "data/subtidal_monitoring/raw/MLPA_kelpforest_site_table.4.csv")) %>%
  janitor::clean_names() %>%
  dplyr::select(site, latitude, longitude, ca_mpa_name_short, mpa_class=site_designation, mpa_designation=site_status, 
                baseline_region)%>%
  #select sites in Carmel and Monterey Bay only
  dplyr::filter(latitude >= 36.46575 & latitude <= 36.64045) %>%
  #drop sites with insufficient data
  dplyr::filter(!(site == "ASILOMAR_DC" |
                    site == "ASILOMAR_UC" |
                    site == "CHINA_ROCK" |
                    site == "CYPRESS_PT_DC" |
                    site == "CYPRESS_PT_UC" |
                    site == "PINNACLES_IN" |
                    site == "PINNACLES_OUT" |
                    site == "PT_JOE" |
                    site == "SPANISH_BAY_DC" |
                    site == "SPANISH_BAY_UC" |
                    site == "BIRD_ROCK"))%>%
  distinct() #remove duplicates

export_dir <- "/Users/jossmith/Desktop/test_data/sst_at_sites"


library(sf)
library(nngeo)

library(sf)

library(sf)


# convert site table to sf object with lon and lat as coordinates
site_sf <- st_as_sf(site_table, coords = c("longitude", "latitude"), crs = 4326)

library(sf)

library(sf)
library(nngeo)

# read in site table
site_table <- read.csv("/Users/jossmith/Desktop/test_data/site_table.csv")

# loop through each year
for (year in 2003:2004) {
  
  # read in .Rds file for current year
  sst_file <- paste0(year, "_sst_daily.Rds")
  sst_data <- readRDS(file.path("/Users/jossmith/Desktop/test_data", sst_file))
  
  # convert sst data to sf object
  sst_sf <- st_as_sf(sst_data, coords = c("longitude", "latitude"), crs = 4326)
  
  # loop through each site in site table
  for (i in 1:nrow(site_table)) {
    
    # create a sf point for the current site
    site_sf <- st_as_sf(site_table[i,], coords = c("longitude", "latitude"), crs = 4326)
    
    # find the closest point in the sst data to the current site
    nearest_sst <- st_nearest_feature(site_sf, sst_sf)
    
    # get the analysed_sst value for the closest point
    if (is.data.frame(nearest_sst)) {
      analysed_sst <- nearest_sst$analysed_sst
    } else {
      analysed_sst <- NA
    }
    
    # save the analysed_sst value to the site table
    site_table[i, paste0("analysed_sst_", year)] <- analysed_sst
    
  }
  
  # save the updated site table for the current year to a new .Rds file
  site_file <- paste0("site_table_", year, ".Rds")
  saveRDS(site_table, file.path("/Users/jossmith/Desktop/test_data/sst_at_sites", site_file))
  
  # clear memory before next iteration
  rm(sst_data, sst_sf)
  
}










