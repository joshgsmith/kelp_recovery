# Read data
################################################################################

# Clear workspace
rm(list = ls())

librarian::shelf(lubridate, tidyverse, sf, janitor, raster)


# Directories
basedir <- "/Volumes/seaotterdb$/kelp_recovery/"
figdir <- here::here("analyses","4patch_drivers","Figures")

#read data

pisco_poly <- st_read(file.path(basedir, "/data/subtidal_monitoring/raw/MONIT_WC_PISCO_SubT_site_polygons/MONIT_WC_PISCO_SubT_site_polygons.shp"))%>%
                janitor::clean_names()%>%
                #filter Monterey region
                filter(region == "CCSR") %>%
                st_transform("+proj=longlat +datum=WGS84")
#plot(pisco_poly)                


vrm_raw <- raster(file.path(basedir, "data/seafloor_data/processed/rug_rast.tif"))

#read state
ca_counties <- st_read(file.path(basedir, "data/gis_data/raw/ca_county_boundaries/s7vc7n.shp")) 

################################################################################
#make it spatial

# Project the raster to CA Teale Albers
vrm_ca <- projectRaster(vrm_raw, crs = "+init=epsg:3310")

vrm_df <- raster::rasterToPoints(vrm_raw)  %>%
              st_as_sf(crs = 3310, coords = c("x","y")) %>%
              st_transform(st_crs = (6414))

################################################################################
#calculate VRM features for each PISCO site

# Reproject pisco_poly to EPSG:3310
pisco_poly <- st_transform(pisco_poly, crs = 3310)

# Reproject vrm_df to EPSG:3310
vrm_df <- st_transform(vrm_df, crs = 3310)


# Join the polygon layer with the raster layer
pisco_raster_join <- sf::st_join(pisco_poly, vrm_df, join = st_intersects)


# Calculate the mean value for each polygon
pisco_mean <- pisco_raster_join %>%
  group_by(site) %>%
  summarise(mean_value = mean(rug_rast, na.rm=TRUE))


