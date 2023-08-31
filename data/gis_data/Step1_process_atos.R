

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
atos_orig <- st_read(file.path(basedir, "gis_data/raw/ATOS/ATOS_mpen_60msw_1000m/ATOS_mpen_60msw_1000m.shp"))

#read state
ca_counties_orig <- st_read(file.path(basedir, "gis_data/raw/ca_county_boundaries/s7vc7n.shp")) 

# Get land
usa <- rnaturalearth::ne_states(country="United States of America", returnclass = "sf")
foreign <- rnaturalearth::ne_countries(country=c("Canada", "Mexico"), returnclass = "sf")

################################################################################
#inspect

st_crs(atos_orig)

plot(atos_orig)

ggplot() +
  geom_sf(data = atos_orig) +
  #add land
  geom_sf(data = ca_counties_orig, fill = "gray", color = "gray80") +
  geom_sf_text(data = atos_orig, aes(label = as.character(ATOS_ID)), 
               color = "black", size = 3, fontface = "bold", nudge_y = 0.001) +
  coord_sf(xlim = c(-121.99, -121.88), ylim = c(36.519, 36.645), crs = 4326) +
  theme_bw()

################################################################################
#create new site name based on upcoast to downcoast area


atos_build1 <- atos_orig %>%
  arrange(ATOS_ID) %>%
  mutate(site_name = row_number()) %>%
  dplyr::select(site_name, geometry)

#check
ggplot() +
  geom_sf(data = atos_build1) +
  #add land
  geom_sf(data = ca_counties_orig, fill = "gray", color = "gray80") +
  geom_sf_text(data = atos_build1, aes(label = as.character(site_name)), 
               color = "black", size = 3, fontface = "bold", nudge_y = 0.001) +
  coord_sf(xlim = c(-121.99, -121.88), ylim = c(36.519, 36.645), crs = 4326) +
  theme_bw()


st_write(atos_build1, file.path(basedir, "gis_data/processed/ATOS/ATOS_mpen_60msw_1000m.shp"))









