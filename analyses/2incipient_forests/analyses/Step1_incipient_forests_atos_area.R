
#Joshua G. Smith; jossmith@mbayaq.org

rm(list=ls())


######
#required packages
librarian::shelf(tidyverse, sf, raster, shiny, tmap)

#set directories 
basedir <- "/Volumes/seaotterdb$/kelp_recovery/data"
figdir <- here::here("analyses","2incipient_forests","figures")
output <- here::here("analyses","2incipient_forests","output")

#read landsat dat
landsat_dat <- st_read(file.path(basedir, "kelp_landsat/processed/monterey_peninsula/landsat_mpen_1984_2023_points_withNAs.shp"))

#read atos area
atos_area <- st_read(file.path(basedir, "gis_data/processed/ATOS/ATOS_mpen_60msw_1000m.shp")) 

#read forage data
forage_dat <- read.csv(file.path(basedir, "foraging_data/processed/foraging_data_2016_2023.csv")) 

#read forage metadata
forage_meta <- read.csv(file.path(basedir, "foraging_data/processed/forage_metadata.csv")) 

################################################################################
#Step 2 - scale biomass to atos area polygon extent 

atos_area <- atos_area %>% st_transform(crs = st_crs(landsat_dat)) 
              
#inspect
# Calculate centroids of polygons
atos_area$centroid <- st_centroid(atos_area)

# Create a data frame for labels and lines
label_data <- data.frame(
  site_name = atos_area$site_name,
  x = st_coordinates(atos_area$centroid)[, "X"],
  y = st_coordinates(atos_area$centroid)[, "Y"],
  xend = st_coordinates(atos_area$centroid)[, "X"],
  yend = st_coordinates(atos_area$centroid)[, "Y"] + 0.002
)

# Scatter plot with jittered site names and lines connecting labels
ggplot(atos_area) +
  geom_sf() +
  geom_segment(
    data = label_data,
    aes(x = x, y = y, xend = xend, yend = yend),
    color = "blue"
  ) +
  ggrepel::geom_label_repel(
    data = label_data,
    aes(x = x, y = y, label = site_name),
    nudge_y = 0.002,
    segment.color = "blue",
    segment.size = 0.5
  ) +
  labs("Longitude", y = "Latitude") +
  theme_minimal()


# Join the point data with the polygon data based on spatial location
joined_data <- st_join(landsat_dat, atos_area)

# Group by year, quarter, and site_name, and summarize the biomass values
summarized_data <- joined_data %>%
  group_by(year, quarter, site_name) %>%
  summarize(total_biomass = sum(biomass, na.rm = TRUE))

################################################################################
#step3 process foraging data

forage_build1 <- forage_dat %>%
  mutate(quarter = ifelse(month == 1 | month == 2 | month == 3, 1,
                          ifelse(month == 4 | month == 5 | month == 6,2,
                                 ifelse(month == 7 | month == 8 | month == 9,3,4))))  %>%
  filter(!(is.na(lat)))%>%
  #make spatial
  st_as_sf(., coords = c("long","lat"), crs=4326)

#join with scan area
forage_build2 <- st_join(forage_build1, atos_area)

#aggregate data by n dives per year and site_name
forage_build3 <- forage_build2 %>%
  group_by(year, site_name) %>%
  summarize(n_dive_success = n()) %>%
  filter(!(is.na(site_name)))

forage_build4 <- st_transform(forage_build3, crs=st_crs(landsat_dat))


#transform to relative effort among sites
annual_total <- forage_build4 %>%
                  group_by(year) %>%
                  dplyr::summarize(total_dives = sum(n_dive_success))

forage_build5 <- forage_build4 %>%
                  st_join(annual_total)%>%
                  mutate(prop_dives = n_dive_success / total_dives)


################################################################################
#prep for plot

#determine baseline average kelp biomass
baseline_average <- summarized_data %>%
  #filter(year >= 2007 & year <= 2013) %>%
  #filter(year >= 2000 & year < 2014) %>% #this sets the baseline period. #use year >= 2016 & year <= 2017
  filter(year < 2014) %>%
  filter(quarter == 3)%>% #filter to quarter 3 for max kelp extent
  group_by(site_name) %>%
  summarize(baseline_avg = mean(total_biomass, na.rm=TRUE),
            baseline_sd = sd(total_biomass, na.rm=TRUE)) %>%
  #drop landsat obs not within a site
  filter(!(is.na(site_name)))


observed_means <- summarized_data %>%
  filter(quarter == 3)%>%
  #filter(site_name %in% c("Point Pinos East", "Point Pinos West", "Pescadero Point West", "Pescadero Point East")) %>%
  group_by(site_name, year) %>%
  summarize(observed_avg = mean(total_biomass, na.rm=TRUE)) 


# Assuming baseline_average and observed_means are data frames
final_data <- observed_means %>%
  #filter(year >= 2000)%>%
  st_join(baseline_average) %>%
  rename(site_name = site_name.x) %>%
  dplyr::select(-site_name.y) %>%
  mutate(deviation = (observed_avg - baseline_avg) / baseline_sd) %>%
  filter(!(is.na(site_name))) %>%
  #clean up
  mutate(deviation = ifelse(deviation %in% c("NaN","Inf"), NA,deviation))



st_write(final_data, file.path(output, "landsat_atos_area_60_1000.geojson"))

