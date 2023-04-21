#Joshua G. Smith; jossmith@mbayaq.org

rm(list=ls())


######
#required packages
librarian::shelf(tidyverse, tidync, sf, raster, terra)

#set directories 
basedir <- "/Volumes/seaotterdb$/kelp_recovery/data"
figdir <- here::here("analyses","figures")

#read landsat dat for MPEN
#landsat_dat <- st_read(file.path(basedir, "kelp_landsat/processed/monterey_peninsula/landsat_mpen_1984_2022_points_withNAs.shp"))

#read state
ca_counties <- st_read(file.path(basedir, "gis_data/raw/ca_county_boundaries/s7vc7n.shp")) 

#read foraging data
#forage_dat <- read.csv(file.path(basedir, "foraging_data/processed/foraging_data_2016_2023.csv")) 
#for geom_spatraster see https://dieghernan.github.io/tidyterra/articles/palettes.html


################################################################################
#step 1 - build base map

# create a bounding box with the desired latitude and longitude values
bbox <- st_bbox(c(xmin = -122.08, ymin = 36.510140, xmax = -121.835531, ymax = 36.670574), crs = st_crs(ca_counties))

# filter the ca_counties object using the bounding box
ca_counties_mpen <- st_intersection(ca_counties, st_as_sfc(bbox))

# check features
ggplot() +
  geom_sf(data = ca_counties_mpen)+
  coord_sf(xlim = c(-121.99, -121.88), ylim = c(36.510140, 36.645))+
  theme_bw()


################################################################################
#step 2 - simulate otter aggregations at one time point

library(sf)

# Define the bounding box
bbox <- st_bbox(c(xmin = -122.08, ymin = 36.510140, xmax = -121.835531, ymax = 36.670574), crs = st_crs(ca_counties))

# filter the ca_counties object using the bounding box
ca_counties_mpen <- st_intersection(ca_counties, st_as_sfc(bbox))

# create a 300 meter buffer around the polygon
buffer <- st_buffer(st_make_valid(ca_counties_mpen), dist = 300)

# get the difference between the buffer and ca_counties_mpen
diff <- st_difference(buffer, ca_counties_mpen)

# create 400 random points inside diff
set.seed(123)
pts <- st_sample(diff, size = 400)



# Plot the objects using ggplot2
ggplot() +
  geom_sf(data = diff) +
  geom_sf(data = pts) +
  coord_sf(xlim = c(-121.99, -121.88), ylim = c(36.519, 36.645)) +
  theme_bw()





