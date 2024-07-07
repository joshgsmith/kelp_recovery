

#Joshua G. Smith; jossmith@mbayaq.org

rm(list=ls())

################################################################################

# This code calculates the annual canopy area for Q3 of each of 123 identified 
# clusters using a rolling average (year t, year t-1, year t+1). The processing
# steps are:

#1: filter the data to quarter 3 for each year
#2: Calculate the total canopy area (sum) for each cluster and each year
#3: Perform the rolling average 
#4: Find the max canopy area for each cluster and each year
#5: Determine the annual perc of max for each year. 

################################################################################

#required packages
librarian::shelf(tidyverse, sf, raster, shiny, tmap, terra, tidyterra, RColorBrewer)

#set directories 
basedir <- "/Volumes/seaotterdb$/kelp_recovery/data"
localdir <- "/Users/jossmith/Documents/Data/landsat"
output <- here::here("analyses","2incipient_forests","output")

#read state
ca_counties_orig <- st_read(file.path(basedir, "gis_data/raw/ca_county_boundaries/s7vc7n.shp")) 

#read landsat raw
landsat_orig <- st_read(file.path(localdir, "/processed/monterey_peninsula/landsat_mpen_1984_2023_points_withNAs.shp"))

#read cluster area
landsat_hclust <- readRDS(file.path(basedir,"/kelp_landsat/processed/landsat_cluster_ID_123.Rds"))

# Get land
usa <- rnaturalearth::ne_states(country="United States of America", returnclass = "sf")
foreign <- rnaturalearth::ne_countries(country=c("Canada", "Mexico"), returnclass = "sf")


################################################################################
# Join landsat data with cluster ID

landsat_build1 <- left_join(landsat_orig, landsat_hclust, by = c("latitude","longitude")) %>%
  #drop sites outside of cluster assigments
  filter(!(is.na(cluster)))


################################################################################
# convert clusters to area

#calculate area of kelp in each cluster
summarized_data <- landsat_build1 %>%
  filter(quarter.x == 3)%>%
  group_by(year.x, cluster)%>%
  summarize(total_area = sum(area.x, na.rm = TRUE))%>%
  #smooth using rolling avg
  arrange(cluster, year.x) %>%
  group_by(cluster) %>%
  mutate(
      roll_area_3 = zoo::rollapply(total_area, width = 3, FUN = mean, 
                                   fill = NA, align = "center", partial = TRUE),
      roll_area_5 = zoo::rollapply(total_area, width = 5, FUN = mean, 
                                   fill = NA, align = "center", partial = TRUE)
    )

#find the max cluster area for the 2004-2013 period
max_cluster_area <- summarized_data %>%
                      filter(year.x > 2004 & year.x < 2014)%>%
                      group_by(cluster)%>%
                      summarize(area_max_3 = max(roll_area_3),
                                area_max_5 = max(roll_area_5))


#join
area_data <- st_join(summarized_data, max_cluster_area) %>%
              select(-cluster.y)%>%
              #clean up
              rename(year = year.x,
                     cluster = cluster.x,
                     )%>%
              mutate(perc_of_max_3 = (roll_area_3 / area_max_3)*100,
                     perc_of_max_5 = (roll_area_5 / area_max_5)*100)%>%
              select(year, cluster, total_area, perc_of_max_3,
                     perc_of_max_5, geometry)
              
plot(area_data %>% filter(year == 2023))
################################################################################
# convert clusters to area



st_write(area_data, file.path(output, "area_data3.geojson")) #last write 16 Feb 2024







