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

# Chlorophyll-a, Aqua MODIS, NPP, 0.0125°, West US, EXPERIMENTAL, 2002-present (MonthlyComposite)

# Build data
################################################################################

# Look up datasets 
datasets_grids <- ed_datasets(which="griddap", url="https://coastwatch.pfeg.noaa.gov/erddap/")
datasets_tables <- ed_datasets(which="tabledap", url="https://coastwatch.pfeg.noaa.gov/erddap/")


# Get data
################################################################################

#set export location
export_path <- file.path(basedir, "data/environmental_data/NPP/processed")

#Set start and end dates
start_date <- as.Date("2003-01-01")
end_date <- as.Date("2020-12-31")

# Download daily NPP
data_info <- info("erdMWchlamday_LonPM180")
data_chunk <- griddap(dataset="erdMWchlamday_LonPM180", 
                      time = c(start_date, end_date),
                      longitude = c(-125, -117), latitude = c(36.500888, 37.122990), #Ano Nuevo to Lobos and offshore 5 Km
                      fields="chlorophyll")

# Combine data frames in data_list
data_list <- list(data_chunk)
data_combined <- do.call(rbind, lapply(data_list, function(x) x$data))

# Remove rows with missing values
data_combined <- data_combined[complete.cases(data_combined), ]

#Make into df and clean
npp_monthly <- data_combined %>%
  mutate(date=gsub("T00:00:00Z", "", time),
         year = lubridate::year(date),
         month = lubridate::month(date),
         day = lubridate::day(date))


# Save data for current year to an Rds file
#saveRDS(npp_monthly, file.path(export_path,"npp_monthly.Rds"))


################################################################################
#calculate NPP at sites

#load monitoring site table
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



###gpt test

library(sf)

site_table_sf <- st_as_sf(site_table, coords = c("longitude", "latitude"), crs = 4326)
npp_monthly_sf <- st_as_sf(npp_monthly, coords = c("longitude", "latitude"), crs = 4326)

# calculate the distance matrix between site_table and npp_monthly
dist_mat <- st_distance(npp_monthly_sf, site_table_sf)

# find the index of the closest point in npp_monthly for each point in site_table
closest_idx <- apply(dist_mat, 1, which.min)

# extract the closest points from npp_monthly
closest_points <- npp_monthly[closest_idx, ]

#add sites
closest_points$site <- site_table$site[closest_idx] 

#find distinct values
filtered_dat <- closest_points %>%
                dplyr::distinct() %>%
                dplyr::rename(npp_lat = `latitude`,
                              npp_long = `longitude`)

#identify sites that didn't join as assign closest neighoboring site
site_index <- left_join(site_table, filtered_dat, by="site") %>%
                dplyr::select(site, npp_long)

#filtered_dat2 <- inner_join(closest_points, site_index, by = c("longitude"="npp_long"))

#apply filtered index to original dataframe
npp_monthly_filtered <- npp_monthly %>% inner_join(site_index, by = c("longitude"="npp_long")) %>% distinct() %>%
  filter(latitude == 36.5)


# Export data
#saveRDS(npp_monthly_filtered, file.path(basedir, "/data/environmental_data/NPP/processed/NPP_at_PISCO_site.Rds"))



