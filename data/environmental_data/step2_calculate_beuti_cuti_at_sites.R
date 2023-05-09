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

# Read CUTI
cuti_orig <- read.csv(file.path(basedir, "/data/environmental_data/CUTI/raw/CUTI_daily.csv"), as.is=T)

#load monitoring site table
site_table <- read.csv(file.path(basedir, "data/subtidal_monitoring/raw/MLPA_kelpforest_site_table.4.csv")) %>%
  janitor::clean_names() %>%
  dplyr::select(site, latitude, longitude, ca_mpa_name_short, mpa_class=site_designation, mpa_designation=site_status, 
                baseline_region)%>%
  distinct() #remove duplicates

# Format CUTI
################################################################################

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
  dplyr::select(year:date, lat_dd_bin, lat_dd_lo, lat_dd_hi, everything())

# Inspect
table(data$lat_dd_bin)

# Export
saveRDS(data, file.path(basedir, "/data/environmental_data/CUTI/processed/1988_2022_cuti_beuti_daily.Rds"))

# 
################################################################################
# Build upwelling time series

sites_bin <- site_table %>% mutate(lat_bin = round(latitude))

site_upwelling <- left_join(data, sites_bin, by=c("lat_dd_bin"="lat_bin")) %>%
                  filter(!(is.na(site)))

# Export data
saveRDS(site_upwelling, file.path(basedir, "/data/environmental_data/CUTI/processed/1988_2022_cuti_beuti_daily_by_PISCO_site.Rds"))
