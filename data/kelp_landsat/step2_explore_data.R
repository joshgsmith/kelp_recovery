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

#read foraging data
forage_dat <- read.csv(file.path(basedir, "foraging_data/processed/foraging_data_2016_2023.csv")) 
#for geom_spatraster see https://dieghernan.github.io/tidyterra/articles/palettes.html

################################################################################
#Step 2 -- process foraging dat for join

forage_build1 <- forage_dat %>%
  mutate(quarter = ifelse(month == 1 | month == 2 | month == 3, 1,
                          ifelse(month == 4 | month == 5 | month == 6,2,
                                 ifelse(month == 7 | month == 8 | month == 9,3,4)))) %>%
  #filter urchin prey only
  filter(prey == "pur") 

focal_patch <-  forage_dat %>%
  mutate(quarter = ifelse(month == 1:3,1,
                          ifelse(month == 4:6,2,
                                 ifelse(month==7:9,3,4)))) %>%
  #filter urchin prey only
  filter(prey == "pur") %>%
  #define focal patch
  group_by(bout) %>%
  dplyr::summarize(n_dives = n()) %>%
  dplyr::mutate(focal_patch = ifelse(n_dives > 3, "yes","no"))


#join

forage_build2 <- left_join(forage_build1, focal_patch, by="bout")



#aggregate data by bout

forage_build3 <- forage_build2 %>%
  group_by(year, quarter, bout, focal_patch) %>%
  summarize(n_prey = sum(preynum),
            lat = mean(lat),
            long = mean(long)
  ) %>%
  filter(!(is.na(lat)))%>%
  #make spatial
  st_as_sf(.,coords = c("long","lat"), crs=4326)



################################################################################
#Test plot to view data structure --- identify coords for MPEN and Carmel

# Get land
#usa <- rnaturalearth::ne_states(country="United States of America", returnclass = "sf")
#foreign <- rnaturalearth::ne_countries(country=c("Canada", "Mexico"), returnclass = "sf")

#Plot Monterey 

plot_dat <- landsat_dat %>% filter ( year == "2019",
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
t_dat <- st_transform(plot_dat, crs=3310) # %>% filter(year == 2019)

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
#  # Plot land
  geom_sf(data=ca_counties, fill="tan",  lwd=0.3) +
  geom_sf(data= forage_build3 %>% filter(year==2019 & quarter == 2) %>% st_transform(crs=3310), 
          aes(color=focal_patch, shape=focal_patch), size=3)+
  scale_color_manual(values=c("blue","red"))+
  # Plot kelp
  labs(x="", y="") +
  # Crop
  coord_sf(xlim = c(-121.93, -121.89), ylim = c(36.60, 36.64), crs=4326) +
  labs(fill = "Biomass \n(kg)")+
  # Theme
  theme_bw() + theme(panel.background = element_rect(fill="#D4FCFF"),
                     panel.grid.major=element_blank())


















