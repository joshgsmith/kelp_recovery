#Joshua G. Smith 
#jossmith@mbayaq.org


# Read data
################################################################################

# Clear workspace
rm(list = ls())

librarian::shelf(rerddap, lubridate, tidyverse, beepr)


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
  # Calculate start and end dates for this year
  chunk_start <- as.Date(paste0(year, "-01-01"))
  chunk_end <- as.Date(paste0(year, "-12-31"))
  
  # Download data for this year
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
  
  sst_daily <- data_combined %>%
    mutate(date=gsub("T00:00:00Z", "", time),
           year = lubridate::year(date),
           month = lubridate::month(date),
           day = lubridate::day(date))
  
  # Create a filename for this year
  filename <- paste0(year, "_sst_daily.Rds")
  
  # Save data for this year to an Rds file
  saveRDS(sst_daily, file.path(export_path,filename))
  
  # Remove data for this year from memory before iterating again
  rm(data_info)
  rm(data_list)
  rm(sst_daily)
  rm(data_combined)
  rm(data_chunk)
}




