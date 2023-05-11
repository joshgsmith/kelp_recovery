#Joshua G. Smith 
#jossmith@mbayaq.org


# Read data
################################################################################

# Clear workspace
rm(list = ls())

librarian::shelf(rerddap, lubridate, tidyverse, beepr, data.table, purrr)


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
#################################################

#set export location
export_path <- file.path("/Users/jossmith/Desktop/test_data")

#Set start and end dates
start_date <- as.Date("2003-01-01")
end_date <- as.Date("2020-12-31")


# Loop over years
for (year in unique(year(start_date):year(end_date))) {
  # Calculate start and end dates for current year
  chunk_start <- as.Date(paste0(year, "-01-01"))
  chunk_end <- as.Date(paste0(year, "-12-31"))
  
  # Download data for current year
  data_info <- info("jplMURSST41")
  data_chunk <- griddap(dataset="jplMURSST41", 
                        time = c(chunk_start, chunk_end),
                        longitude = c(-125, -117), latitude = c(36.500888, 37.122990), #Ano Nuevo to Lobos
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
#################################################

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
#################################################

sst_dat <- readRDS("2003_2022_daily_SST.Rds")

names(sst_dat)





