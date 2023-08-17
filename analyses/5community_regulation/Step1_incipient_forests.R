
#Joshua G. Smith; jossmith@mbayaq.org

rm(list=ls())


######
#required packages
librarian::shelf(tidyverse, sf, raster, shiny, tmap)

#set directories 
basedir <- "/Volumes/seaotterdb$/kelp_recovery/data"
figdir <- here::here("analyses","figures")

#read landsat dat
landsat_dat <- st_read(file.path(basedir, "kelp_landsat/processed/monterey_peninsula/landsat_mpen_1984_2022_points_withNAs.shp"))

#read state
ca_counties <- st_read(file.path(basedir, "gis_data/raw/ca_county_boundaries/s7vc7n.shp")) 

#read otter scan area
scan_area <- st_read(file.path(basedir, "gis_data/raw/otter_scan_area/TrackingMaps_POLY.shp")) 

################################################################################
#Step 1 - build stacked raster

rast_raw <- landsat_dat 

# transform landsat data to Teale Albers
rast_build1 <- st_transform(rast_raw, crs = 3310) %>% 
  mutate(biomass = ifelse(biomass == 0, NA, biomass))  %>% filter (year >= 2007) #reduce years if memory issue

#define blank raster
r <- terra::rast(rast_build1, res=30)

# Get unique quarters and years
quarters <- unique(rast_build1$quarter)
years <- unique(rast_build1$year)

# Loop through each quarter and year
raster_list <- list()
for (y in years) {
  for (q in quarters) {
    # Subset rast_build1 for current quarter and year
    subset <- rast_build1 %>% filter(year == y, quarter == q)
    # Rasterize subset
    rasterized <- rasterize(subset, r, field = "biomass", fun = mean)
    # Add year and quarter as part of the SpatRaster object name
    name <- paste0("year_", y, "_quarter_", q)
    names(rasterized) <- name
    # Add raster to list
    raster_list <- append(raster_list, rasterized)
  }
}

# Stack rasters in list
stacked_raster <- stack(raster_list)

#check names
names(stacked_raster)

str(stacked_raster)

################################################################################
#Step 2 - scale biomass to scan area polygon extent 

#glace at scan area
plot(scan_area)

scan_area <- scan_area %>% rename(site_name = Name) %>%  st_transform(crs = crs(stacked_raster))

# Get the bounding box of the transformed polygons
polygon_bbox <- st_bbox(scan_area)

# Set the extent of the raster to match the bounding box of polygons
extent(stacked_raster) <- polygon_bbox

# Initialize an empty list to store the sum values
sum_values <- list()

# Initialize an empty list to store the sum matrices
sum_matrices <- list()

# Loop through each layer in the raster stack
for (i in 1:nlayers(stacked_raster)) {
  layer <- stacked_raster[[i]]
  
  # Initialize an empty matrix to store the sums for each polygon
  polygon_sums <- matrix(0, nrow(scan_area), 1)
  
  # Loop through each polygon in scan_area
  for (j in 1:nrow(scan_area)) {
    polygon <- scan_area[j, ]
    
    # Crop the layer to the polygon extent
    cropped_layer <- crop(layer, polygon)
    
    # Calculate the sum of values within the polygon extent
    polygon_sum <- sum(values(cropped_layer), na.rm = TRUE)
    
    polygon_sums[j, 1] <- polygon_sum
  }
  
  # Add the polygon sums matrix for the current layer to the sum_matrices list
  sum_matrices[[i]] <- polygon_sums
}

# Construct the final sum matrix by binding the matrices horizontally
final_sum_matrix <- do.call(cbind, sum_matrices)

# Create the final data frame by combining site names and sum matrix
sum_df <- data.frame(scan_area$site_name, final_sum_matrix)

# Set column names
colnames(sum_df) <- c("site_name", names(stacked_raster))


################################################################################
#Step 3 - transpose data

kelp_area_build1 <- sum_df %>%
  pivot_longer(
    cols = 2:ncol(.),
    names_to = "Layer",
    values_to = "kelp_biomass"
  ) %>%
  #bust apart the meta data
  mutate(
    year = sub("year_(\\d+)_quarter_\\d+", "\\1", Layer),
    quarter = sub("year_\\d+_quarter_(\\d+)", "\\1", Layer)
  ) %>%
  dplyr::select(-Layer)
  

################################################################################
#plot

library(dplyr)
library(ggplot2)


baseline_average <- kelp_area_build1 %>%
  #filter(year >= 2007 & year <= 2013) %>%
  filter(year >= 2015 & year <= 2017) %>%
  group_by(site_name) %>%
  summarize(baseline_avg = mean(kelp_biomass))

observed_means <- kelp_area_build1 %>%
  filter(site_name %in% c("Point Pinos East", "Point Pinos West", "Pescadero Point West", "Pescadero Point East")) %>%
  group_by(site_name, year) %>%
  summarize(observed_avg = mean(kelp_biomass))

final_data <- observed_means %>%
  left_join(baseline_average, by = "site_name") %>%
  mutate(deviation = (observed_avg - baseline_avg) / sd(observed_avg))

# Plotting
ggplot(final_data, aes(x = year, y = deviation, group = site_name)) +
  geom_line() +
  geom_point() +
  facet_wrap(~site_name, ncol = 1, scales = "free_y") +
  labs(x = "Year", y = "Deviation (Standard Deviations from Baseline)", title = "Observed Annual Mean vs. Baseline Deviation by Site") +
  theme_minimal()

####to do: try setting the reference (baseline) period as 2015-2017ish so
#that we can better track recovery
