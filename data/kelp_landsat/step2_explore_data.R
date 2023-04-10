#Joshua G. Smith, jossmith@mbayaq.org

rm(list=ls())

#required packages
librarian::shelf(tidyverse, tidync, sf, rnaturalearth, usethis, raster)

#set directories 
basedir <- "/Volumes/seaotterdb$/kelp_recovery/data"
figdir <- here::here("analyses","figures")

#read landsat dat
landsat_dat <- st_read(file.path(basedir, "kelp_landsat/processed/monterey_peninsula/landsat_mpen_1984_2022_points.shp"))

#read state
ca_counties <- st_read(file.path(basedir, "gis_data/raw/ca_county_boundaries/s7vc7n.shp")) 


#for geom_spatraster see https://dieghernan.github.io/tidyterra/articles/palettes.html

################################################################################
#Test plot to view data structure --- identify coords for MPEN and Carmel

# Get land
#usa <- rnaturalearth::ne_states(country="United States of America", returnclass = "sf")
#foreign <- rnaturalearth::ne_countries(country=c("Canada", "Mexico"), returnclass = "sf")

#Plot Monterey 

plot_dat <- landsat_dat %>% filter (quarter == "2" &
                                          year == "2019" &
                                          latitude >= 36.510140 &
                                          latitude <= 36.670574) %>%
  filter(!(biomass == "0" | is.na(biomass)))


##Entire MPEN coords

mpen <- ggplot() +
  # Plot land
  geom_sf(data=ca_counties, fill="tan",  lwd=0.3) +
  # Plot kelp
  geom_sf(data=plot_dat, aes(color = biomass), size = 1) +
  # Labels
  labs(x="", y="") +
  # Crop
  coord_sf(xlim = c(-122, -121.87), ylim = c(36.50, 36.65)) +
  # Theme
  theme_bw()

mpen

##PG coords
pg <- ggplot() +
  # Plot land
  geom_sf(data=ca_counties, fill="tan",  lwd=0.3) +
  # Plot kelp
  geom_sf(data=plot_dat, aes(color = biomass), size = 1) +
  # Labels
  labs(x="", y="") +
  # Crop
  coord_sf(xlim = c(-121.93, -121.89), ylim = c(36.60, 36.64)) +
  # Theme
  theme_bw()

pg


##Cypress coast coords
cc <- ggplot() +
  # Plot land
  geom_sf(data=ca_counties, fill="tan",  lwd=0.3) +
  # Plot kelp
  geom_sf(data=plot_dat, aes(color = biomass), size = 1) +
  # Labels
  labs(x="", y="") +
  # Crop
  coord_sf(xlim = c(-121.99, -121.93), ylim = c(36.56, 36.64)) +
  # Theme
  theme_bw()


cc


##Carmel Bay
cb <- ggplot() +
  # Plot land
  geom_sf(data=ca_counties, fill="tan",  lwd=0.3) +
  # Plot kelp
  geom_sf(data=plot_dat, aes(color = biomass), size = 1) +
  # Labels
  labs(x="", y="") +
  # Crop
  coord_sf(xlim = c(-121.96, -121.92), ylim = c(36.515, 36.569)) +
  # Theme
  theme_bw()

cb


################################################################################
#Rasterize

#transform landsat data to Teale Albers
t_dat <- st_transform(plot_dat, crs=3310)

#create grid
r <- terra::rast(t_dat, res=30)  # Builds a blank raster of given dimensions and resolution  
vr <- rasterize(t_dat, r,"biomass", resolution = 30) 


#examine output

#MPEN
ggplot() +
  tidyterra::geom_spatraster(data=vr, na.rm=T) +
  tidyterra::scale_fill_whitebox_c(
    palette = "viridi",
    labels = scales::label_number(#suffix = "Density",
    )
  )+
  # Plot land
  geom_sf(data=ca_counties1, fill="tan",  lwd=0.3) +
  # Plot kelp
  labs(x="", y="") +
  # Crop
  coord_sf(xlim = c(-121.96, -121.92), ylim = c(36.515, 36.569), crs=4326) +
  labs(fill = "Biomass \n(kg)")+
  # Theme
  theme_bw() + theme(panel.background = element_rect(fill="#D4FCFF"),
                     panel.grid.major=element_blank())


















