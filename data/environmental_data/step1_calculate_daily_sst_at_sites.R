
# Read data
################################################################################

# Clear workspace
rm(list = ls())

librarian::shelf(rerddap, lubridate, tidyverse)


# Directories
basedir <- "/Volumes/seaotterdb$/kelp_recovery/"
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

# Get data 

# Set start and end dates
start_date <- as.Date("2002-06-16")
end_date <- as.Date("2022-12-31")

# Calculate number of days between start and end dates
days_diff <- as.numeric(end_date - start_date)

# Define chunk size 
chunk_size <- days_diff / 20

# Calculate number of chunks
num_chunks <- ceiling(days_diff / chunk_size)

# Create empty list to store data
data_list <- list()

# Loop over chunks
for (i in seq(0, days_diff, by = chunk_size)) {
  # Calculate start and end dates for this chunk
  chunk_start <- start_date + i
  chunk_end <- min(chunk_start + chunk_size - 1, end_date)
  
  # Download data for this chunk
  data_info <- info("jplMURSST41")
  data_chunk <- griddap(dataset="jplMURSST41", 
                        time = c(chunk_start, chunk_end),
                        longitude = c(-125, -117), latitude = c(34.48084, 38.95117), #vandenberg to Pt. Arena
                        fields="analysed_sst")
  
  # Append data to list
  data_list[[i/chunk_size + 1]] <- data_chunk
}


# Combine data frames in data_list
data_combined <- do.call(rbind, lapply(data_list, function(x) x$data))

# Remove rows with missing values
data_combined <- data_combined[complete.cases(data_combined), ]


sst_daily <- data_combined %>%
                mutate(date=gsub("T00:00:00Z", "", time),
                       year = lubridate::year(date),
                       month = lubridate::month(date),
                       day = lubridate::day(date))


# Inspect
head(sst_daily)
str(sst_daily)


# Export
#saveRDS(data_df1, file=file.path(outputdir, "2002_2022_mursst_monthly.Rds"))


# Convert to raster
#################################################

# Test dataframe
test <- data_df1 %>% 
  filter(date=="2002-06-16") %>% 
  mutate(lat_dd=round(lat_dd, 2),
         long_dd=round(long_dd, 2))

# Convert df to raster brick
data_ras <- freeR::df2brick(data_df1, x_col = "long_dd", y_col="lat_dd", z_col = "sst_c", layer_col = "date")

# Project
data_ras_prj <- data_ras 
raster::crs(data_ras_prj) <- "+proj=longlat +datum=WGS84"

# Chekc
raster::plot(data_ras_prj, 3)

# Export
raster::writeRaster(data_ras_prj, file.path(outputdir, "2002_2022_mursst_monthly_raster.grd"))






