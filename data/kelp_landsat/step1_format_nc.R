#Joshua G. Smith, jossmith@mbayaq.org
# nc processing adapted from Julien Brun, brun@nceas.ucsb.edu

# Data Source on EDI
# https://portal.edirepository.org/nis/mapbrowse?packageid=knb-lter-sbc.74.17


rm(list=ls())

#set memory #### USE R_MAX_VSIZE=100Gb; 
## BE SURE TO SAVE .Renviron
usethis::edit_r_environ()

#restart R 
.rs.restartR()

# install.packages("librarian")
librarian::shelf(tidyverse, tidync, sf, rnaturalearth, usethis, raster)

#set directories 
datin <- "/Users/jossmith/Documents/Data/landsat"
#datin <- "/Volumes/seaotterdb$/kelp_recovery/data/kelp_landsat/raw"
datout <- "/Volumes/seaotterdb$/kelp_recovery/data/kelp_landsat/processed"
figdir <- here::here("analyses","figures")

landsat_dat <- "LandsatKelpBiomass_2023_Q2_withmetadata.nc"

# read the data in
landsat_raw <- tidync(file.path(datin, landsat_dat)) 

# Transform the biomass grid into a data frame
landsat_df <- landsat_raw %>% 
  activate() %>%
  hyper_tibble(force = TRUE) 

# Transform the time, year, qtr grid into a data frame
## Note: time is unique identifier
kelp_year <- landsat_raw %>%
  activate("D1") %>%
  hyper_tibble() 


# Transform the latitude grid into a data frame
kelp_lat <- landsat_raw %>%
  activate("latitude") %>%
  hyper_tibble()

# Transform the longitude grid into a data frame
kelp_lon <- landsat_raw %>%
  activate("longitude") %>%
  hyper_tibble()

# Join lat and lon
kelp_latlon <- left_join(kelp_lat, kelp_lon, by = "station") %>%
  relocate(station, .before=everything())

# Join landsat df with kelp latlong
kelp_df_latlon <- left_join(landsat_df, kelp_latlon, by = "station") %>%
  relocate(station, .before=everything())

#check join
head(kelp_df_latlon)

# Join df with time
kelp_processed <- left_join(kelp_df_latlon, kelp_year, by = "time") %>%
  relocate(station, time, year, quarter, .before=everything())

#check join
head(kelp_processed)
nrow(kelp_processed)

# transform it to sf object
kelp_landsat_sf <- kelp_processed %>%
  st_as_sf(coords=c("longitude", "latitude"), crs=4326, remove=FALSE)

rm(kelp_processed)
rm(kelp_df_latlon)
rm(kelp_latlon)
rm(kelp_lon)
rm(kelp_lat)
rm(landsat_df)

################################################################################
#inspect
kelp_2008 <- kelp_landsat_sf %>% filter(year == 2008 & quarter == 3) %>%
                              mutate(biomass = ifelse(biomass > 0,1,0))

#MPEN
ggplot() +
  # Plot land
  #geom_sf(data=foreign, fill="grey90", color="white", lwd=0.3) +
 # geom_sf(data=usa, fill="grey90", color="white", lwd=0.3) +
  # Plot kelp
  geom_sf(data=kelp_2008, aes(color = area), size = 1) +
  # Labels
  labs(x="", y="") +
  # Crop
  coord_sf(xlim = c(-122.08, -121.835531), ylim = c(36.510140, 36.670574)) +
  # Theme
  theme_bw()

#GFNMS
ggplot() +
  # Plot land
  #geom_sf(data=foreign, fill="grey90", color="white", lwd=0.3) +
  #geom_sf(data=usa, fill="grey90", color="white", lwd=0.3) +
  # Plot kelp
  geom_sf(data=kelp_2008, aes(color = area), size = 1) +
  # Labels
  labs(x="", y="") +
  # Crop
  coord_sf(xlim = c(-123.753617, -123.238646), ylim = c(38.508139, 38.953937))+
  # Theme
  theme_bw()


################################################################################
#export processed data by focal regions

#Monterey peninsula
kelp_landsat_export_MPEN <- kelp_landsat_sf %>%
                      #filter Monterey Peninsula
                      filter(
                        #Malpaso creek = southernmost point
                        latitude >= 36.481281 &
                        #Marina SB = northernmost point
                        latitude <= 36.696714)
                        
#Greater Farallones
kelp_landsat_export_GF <- kelp_landsat_sf %>%
  #filter Monterey Peninsula
  filter(
    #Ft. Ross = southernmost point
    latitude >= 38.508139 &
      #Pt. Arena = northernmost point
      latitude <= 38.953937)

#Trinidad
kelp_landsat_export_trin <- kelp_landsat_sf %>%
  #filter Monterey Peninsula
  filter(
    #Ft. Ross = southernmost point
    latitude >= 40.732094 &
      #Pt. Arena = northernmost point
      latitude <= 41.326140)


################################################################################
#Test plot to view data structure

# Get land
usa <- rnaturalearth::ne_states(country="United States of America", returnclass = "sf")
foreign <- rnaturalearth::ne_countries(country=c("Canada", "Mexico"), returnclass = "sf")

#--------------------------Plot Monterey--------------------------------------

plot_mpen <- kelp_landsat_export_MPEN %>% filter (quarter == "3" &
                           year == "2022")

mpen <- ggplot() +
  # Plot land
  #geom_sf(data=foreign, fill="grey90", color="white", lwd=0.3) +
  geom_sf(data=usa, fill="grey90", color="white", lwd=0.3) +
  # Plot kelp
  geom_sf(data=plot_mpen, aes(color = area), size = 1) +
  # Labels
  labs(x="", y="") +
  # Crop
  coord_sf(xlim = c(-122.08, -121.835531), ylim = c(36.510140, 36.670574)) +
  # Theme
  theme_bw()

mpen

#--------------------------Plot Monterey--------------------------------------

plot_GF <- kelp_landsat_export_GF %>% filter (quarter == "3" &
                                                    year == "2008")

GF <- ggplot() +
  # Plot land
  #geom_sf(data=foreign, fill="grey90", color="white", lwd=0.3) +
  geom_sf(data=usa, fill="grey90", color="white", lwd=0.3) +
  # Plot kelp
  geom_sf(data=plot_GF, aes(color = biomass), size = 1) +
  # Labels
  labs(x="", y="") +
  # Crop
  coord_sf(xlim = c(-123.753617, -123.238646), ylim = c(38.508139, 38.953937))+
  # Theme
  theme_bw()

GF



################################################################################
#test plot

ggplot() +
  # Plot land
  #geom_sf(data=foreign, fill="grey90", color="white", lwd=0.3) +
 # geom_sf(data=usa, fill="grey90", color="white", lwd=0.3) +
  # Plot kelp
  geom_sf(data=plot_2008, aes(color = biomass), size = 1) +
  # Labels
  labs(x="", y="") +
  # Crop
  #coord_sf(xlim = c(-123.753617, -123.238646), ylim = c(38.508139, 38.953937))+
  # Theme
  theme_bw()

################################################################################
#Export
st_write(kelp_landsat_export_MPEN, file.path(datin, "processed/monterey_peninsula/landsat_mpen_1984_2023_points_withNAs.shp"))
st_write(kelp_landsat_export_GF, file.path(datin, "processed/gfnms_sonoma/landsat_mpen_1984_2023_points_withNAs.shp"))
st_write(kelp_landsat_export_trin, file.path(datin, "processed/trinidad/landsat_mpen_1984_2023_points_withNAs.shp"))


