

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
  mutate(site_num = row_number()) %>%
  dplyr::select(site_num, geometry) %>%
  #change site name to something more useful
  mutate(site_name = case_when(
    site_num == 1 ~ "Coast Guard Pier",
    site_num == 2 ~ "Cannery Row",
    site_num == 3 ~ "Hopkins",
    site_num == 4 ~ "Lovers Point",
    site_num == 5 ~ "Otter Point East",
    site_num == 6 ~ "Piños East",
    site_num == 7 ~ "Piños",
    site_num == 8 ~ "Piños West",
    site_num == 9 ~ "Asilomar",
    site_num == 10 ~ "Spanish Bay",
    site_num == 11 ~ "Moss Beach",
    site_num == 12 ~ "Point Joe East",
    site_num == 13 ~ "Point Joe West",
    site_num == 14 ~ "China Rock",
    site_num == 15 ~ "Bird Rock North",
    site_num == 16 ~ "Bird Rock South",
    site_num == 17 ~ "Fanshell Beach",
    site_num == 18 ~ "Cypress Point Golf Course",
    site_num == 19 ~ "Sunset Point",
    site_num == 20 ~ "Lone Cypress",
    site_num == 21 ~ "Pescadero East",
    site_num == 22 ~ "Pescadero West",
    site_num == 23 ~ "Stillwater Cove",
    site_num == 24 ~ "Carmel Beach North",
    site_num == 25 ~ "Scenic Road",
    site_num == 26 ~ "Butterfly House",
    site_num == 27 ~ "Carmel River Beach",
    site_num == 28 ~ "Monastery",
    site_num == 29 ~ "Mono-Lobo",
    site_num == 30 ~ "Whaler's Cove",
    site_num == 31 ~ "Sea Lion Point",
    site_num == 32 ~ "South Lobos",
    TRUE ~ NA
  ))





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









