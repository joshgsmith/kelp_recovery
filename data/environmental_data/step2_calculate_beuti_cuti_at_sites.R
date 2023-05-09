##Joshua G. Smith
#May 8, 2023

rm(list=ls())
librarian::shelf(tidyverse)

#processing script adapted from: https://github.com/NCEAS/ca-mpa/blob/main/data/environmental/Step1_format_beuti_cuti_data.R

################################################################################


# Set Directories
basedir <-"/Volumes/seaotterdb$/kelp_recovery/"
figdir <- here::here("analyses","4patch_drivers","Figures")

# Downloaded from: https://mjacox.com/upwelling-indices/

# Format CUTI
################################################################################

# Read CUTI
cuti_orig <- read.csv(file.path(basedir, "/data/environmental_data/CUTI/raw/CUTI_daily.csv"), as.is=T)

# Format CUTI
cuti <- cuti_orig %>% 
  # Gather
  gather(key="lat_dd_bin", value="cuti", 4:ncol(.)) %>% 
  # Format lat bin
  mutate(lat_dd_bin=gsub("X|N", "", lat_dd_bin) %>% as.numeric(.)) %>% 
  # Build date
  mutate(date=paste(year, month, day, sep="-") %>% lubridate::ymd()) %>% 
  # Arrange
  dplyr::select(year, month, day, date, lat_dd_bin, cuti)


# Format BEUTI
################################################################################

# Read CUTI
beuti_orig <- read.csv(file.path(basedir, "/data/environmental_data/BEUTI/raw/BEUTI_daily.csv"), as.is=T)

# Format CUTI
beuti <- beuti_orig %>% 
  # Gather
  gather(key="lat_dd_bin", value="beuti", 4:ncol(.)) %>% 
  # Format lat bin
  mutate(lat_dd_bin=gsub("X|N", "", lat_dd_bin) %>% as.numeric(.)) %>% 
  # Build date
  mutate(date=paste(year, month, day, sep="-") %>% lubridate::ymd()) %>% 
  # Arrange
  dplyr::select(year, month, day, date, lat_dd_bin, beuti)


# Merge data
################################################################################

# Build data
data <- cuti %>% 
  left_join(beuti) %>% 
  # Add lat lo/hi
  mutate(lat_dd_lo=lat_dd_bin-0.5,
         lat_dd_hi=lat_dd_bin+0.5) %>% 
  # Arrange
  dplyr::select(year:date, lat_dd_bin, lat_dd_lo, lat_dd_hi)

# Inspect
table(data$lat_dd_bin)

# Export
saveRDS(data, file.path(basedir, "/data/environmental_data/CUTI/processed/1988_2022_cuti_beuti_daily.Rds"))

# Merge data
################################################################################

