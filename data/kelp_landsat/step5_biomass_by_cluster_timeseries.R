

#Joshua G. Smith; jossmith@mbayaq.org

rm(list=ls())

######
######
#required packages
librarian::shelf(tidyverse, sf, raster, shiny, tmap, terra, tidyterra, RColorBrewer)

#set directories 
basedir <- "/Volumes/seaotterdb$/kelp_recovery/data"
output <- here::here("analyses","2incipient_forests","output")

#read state
ca_counties_orig <- st_read(file.path(basedir, "gis_data/raw/ca_county_boundaries/s7vc7n.shp")) 

#read landsat raw
landsat_orig <- st_read(file.path(basedir, "kelp_landsat/processed/monterey_peninsula/landsat_mpen_1984_2023_points_withNAs.shp"))

#read cluster area
landsat_hclust <- readRDS(file.path(basedir,"/kelp_landsat/processed/landsat_cluster_ID.Rds"))

# Get land
usa <- rnaturalearth::ne_states(country="United States of America", returnclass = "sf")
foreign <- rnaturalearth::ne_countries(country=c("Canada", "Mexico"), returnclass = "sf")

################################################################################
# Join landsat data with cluster ID

landsat_build1 <- left_join(landsat_orig, landsat_hclust, by = c("latitude","longitude")) %>%
                    #drop sites outside of cluster assigments
                    filter(!(is.na(site_name)))


################################################################################
# summarize

# Group by year, quarter, and site_name, and summarize the biomass values
# We want the total annual bioamss for each cluster
summarized_data <- landsat_build1 %>%
  #filter to year 1990 and beyond for Q3 only
  filter(year >= 1990 & quarter == 3)%>%
  group_by(year, quarter, site_num, site_name) %>%
  summarize(total_biomass = sum(biomass, na.rm = TRUE))

#determine baseline average kelp biomass
baseline_average <- summarized_data %>%
  filter(year < 2014) %>%
  filter(quarter == 3)%>% #filter to quarter 3 
  group_by(site_num, site_name) %>%
  summarize(baseline_avg = mean(total_biomass, na.rm=TRUE),
            baseline_sd = sd(total_biomass, na.rm=TRUE)) 

#inspect
#View(baseline_average)

# calculate departures in standard deviation from baseline avg
final_data <- summarized_data %>%
  st_join(baseline_average) %>%
  rename(site_num = site_num.x,
         site_name = site_name.x) %>%
  dplyr::select(-site_name.y, -site_num.y) %>%
  mutate(deviation = (total_biomass - baseline_avg) / baseline_sd)

st_write(final_data, file.path(output, "landsat_cluster_area.geojson"))







