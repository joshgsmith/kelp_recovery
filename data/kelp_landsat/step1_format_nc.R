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
datin <- "/Volumes/seaotterdb$/kelp_recovery/data/kelp_landsat/raw"
datout <- "/Volumes/seaotterdb$/kelp_recovery/data/kelp_landsat/processed"
figdir <- here::here("analyses","figures")

landsat_dat <- "LandsatKelpBiomass_2023_Q1_withmetadata.nc"

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

################################################################################
#export processed data

kelp_landsat_export <- kelp_landsat_sf %>%
                      #drop NAs and 0s to shorten file size (include if historic cover is of interest)
                      #filter(!(is.na(biomass) | biomass == "0")) %>%
                      #filter Monterey Peninsula
                      filter(
                        #Malpaso creek = southernmost point
                        latitude >= 36.481281 &
                        #Marina SB = northernmost point
                        latitude <= 36.696714)
                        
#nrow(kelp_landsat_export)

st_write(kelp_landsat_export, file.path(datout, "monterey_peninsula/landsat_mpen_1984_2023_points_withNAs.shp"))


################################################################################
#Test plot to view data structure

# Get land
usa <- rnaturalearth::ne_states(country="United States of America", returnclass = "sf")
foreign <- rnaturalearth::ne_countries(country=c("Canada", "Mexico"), returnclass = "sf")

#Plot Monterey 

plot_dat <- kelp_landsat_export %>% filter (quarter == "2" &
                           year == "2019" &
                             latitude >= 36.510140 &
                             latitude <= 36.670574) # %>%
                 # filter(!(biomass == "0" | is.na(biomass)))
                    

mpen <- ggplot() +
  # Plot land
  geom_sf(data=foreign, fill="grey90", color="white", lwd=0.3) +
  geom_sf(data=usa, fill="grey90", color="white", lwd=0.3) +
  # Plot kelp
  geom_sf(data=plot_dat, aes(color = biomass), size = 1) +
  # Labels
  labs(x="", y="") +
  # Crop
  coord_sf(xlim = c(-122.08, -121.835531), ylim = c(36.510140, 36.670574)) +
  # Theme
  theme_bw()

mpen




