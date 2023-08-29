

#Joshua G. Smith; jossmith@mbayaq.org

rm(list=ls())


######
#required packages
librarian::shelf(tidyverse, sf, raster)

#set directories 
basedir <- "/Volumes/seaotterdb$/kelp_recovery/data"
figdir <- here::here("analyses","2incipient_forests","figures")
output <- here::here("analyses","2incipient_forests","output")

#read USGS Atos
atos_orig <- st_read(file.path(basedir, "gis_data/raw/ATOS/ATOS_polygon_teale83/ATOS_polygon_teale83.shp"))

#read state
ca_counties_orig <- st_read(file.path(basedir, "gis_data/raw/ca_county_boundaries/s7vc7n.shp")) 

# Get land
usa <- rnaturalearth::ne_states(country="United States of America", returnclass = "sf")
foreign <- rnaturalearth::ne_countries(country=c("Canada", "Mexico"), returnclass = "sf")

################################################################################
#inspect

st_crs(atos_orig)

################################################################################
#select Monterey Bay as focal region

atos_mon <- atos_orig %>%
  st_transform(st_crs(ca_counties_orig)) %>%
  # Crop to Monterey region
  st_crop(st_bbox(c(ymin = 36.519, xmin = -121.99, ymax = 36.7, xmax = -121.88)))

#inspect
plot(atos_mon)

ggplot() +
  geom_sf(data = atos_mon) +
  #add land
  geom_sf(data = ca_counties_orig, fill = "gray", color = "gray80") +
  coord_sf(xlim = c(-121.99, -121.88), ylim = c(36.519, 36.645), crs = 4326) +
  theme_bw()





