#Joshua G. Smith; jossmith@mbayaq.org

rm(list=ls())


######
#required packages
librarian::shelf(tidyverse, tidync, sf, raster, terra)

#set directories 
basedir <- "/Volumes/seaotterdb$/kelp_recovery/data"
figdir <- here::here("analyses","figures")

#read landsat dat
landsat_dat <- st_read(file.path(basedir, "kelp_landsat/processed/monterey_peninsula/landsat_mpen_1984_2022_points_withNAs.shp"))

#read state
ca_counties <- st_read(file.path(basedir, "gis_data/raw/ca_county_boundaries/s7vc7n.shp")) 

#read foraging data
forage_dat <- read.csv(file.path(basedir, "foraging_data/processed/foraging_data_2016_2023.csv")) 
#for geom_spatraster see https://dieghernan.github.io/tidyterra/articles/palettes.html

################################################################################
#Step 1 - rasterize by year and quarter

rast_raw <- landsat_dat

# transform landsat data to Teale Albers
rast_build1 <- st_transform(rast_raw, crs = 3310) %>% 
  mutate(biomass = ifelse(biomass == 0, NA, biomass)) # %>% filter (year == 2016 | year == 2017 | year == 2018)

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

# test plot random layer in the stack
plot(stacked_raster[[75]])


################################################################################
#Step 2 - process historical extent as layer

###code below is deprecated

#define 0 cover as historical kelp footprint
na_dat <- landsat_dat %>% filter (latitude >= 36.510140 &
                                    latitude <= 36.670574, 
                                  biomass == 0)

kelp_historic <- st_transform(na_dat, crs=3310) # %>% mutate(biomass = ifelse(biomass==0,1,NA))


kelp_na <- rasterize(kelp_historic, r,"biomass", resolution = 30) 




