
#Joshua G. Smith; jossmith@mbayaq.org

rm(list=ls())


######
#required packages
librarian::shelf(tidyverse, sf, zoo)

#set directories 
basedir <- "/Volumes/seaotterdb$/kelp_recovery/data"
gisdir <- "/Volumes/seaotterdb$/kelp_recovery/data/gis_data"

#read long term data
lt_orig <- readxl::read_xlsx(file.path(basedir,"intertidal_monitoring/raw/long_term_surveys/Phototrans_summary.xlsx"),sheet=1) %>%
            janitor::clean_names()

#read species table
spp_tab <- readxl::read_xlsx(file.path(basedir,"intertidal_monitoring/raw/spp_codes.xlsx"),sheet=1) %>%
           janitor::clean_names()

#read state
ca_counties_orig <- st_read(file.path(basedir, "gis_data/raw/ca_county_boundaries/s7vc7n.shp")) 


################################################################################
#step 1 - join taxonomy and select vars

lt_build1 <- lt_orig %>% left_join(spp_tab, by = c("lumping_code"="species_code"))

#inspect
unique(lt_build1$marine_species_name)
names(lt_build1)

#select vars

lt_build2 <- lt_build1 %>% dplyr::select(marine_site_name, site_code, latitude,
                                         longitude, marine_common_year, sampling_period,
                                         season_name, target_assemblage, lumping_code,
                                         marine_species_name, 
                                         average_percent_cover, num_plots_sampled,
                                         stddev, stderr) %>% 
              #select california
              filter(latitude < 42)

################################################################################
#step 2 - inspect distribution of sites

marine_sites <- lt_build2 %>% select(latitude, longitude) %>% st_as_sf(coords = c("longitude", "latitude"), crs=4326)



# Plot marine_sites
marine_plot <- ggplot() +
  geom_sf(data = marine_sites, color = "blue", size = 1)+
  geom_sf(data = ca_counties_orig, fill = "gray80", size = 0.5)+
  #pismo to santa cruz
  coord_sf(ylim = c(35.084635, 37.124608), crs = 4326)
  

# Print the combined plot
print(marine_plot)


# Plot marine_sites MPEN
marine_plot <- ggplot() +
  geom_sf(data = marine_sites, color = "blue", size = 4)+
  geom_sf(data = ca_counties_orig, fill = "gray80", size = 0.5)+
  #pismo to santa cruz
  coord_sf(xlim = c(-121.99, -121.88), ylim = c(36.45, 36.645), crs = 4326)


# Print the combined plot
print(marine_plot)


################################################################################
#step 3 - gather data for MPEN


lt_build3 <- lt_build2 %>% dplyr::filter(latitude >= 36.45 & latitude <= 36.645) %>%
              #select focal species
              filter(marine_species_name %in% c("mytilus californianus",
                                                "pisaster ochraceus"))





# Calculate the mean percent cover for each species by year
mean_percent_cover <- lt_build3 %>%
  group_by(marine_common_year, marine_species_name) %>%
  summarize(mean_percent_cover = mean(average_percent_cover, na.rm = TRUE))

# Create a ggplot
ggplot(mean_percent_cover, aes(x = marine_common_year, y = mean_percent_cover, color = marine_species_name)) +
  geom_line() +
  geom_point() +
  labs(x = "Marine Common Year", y = "Mean Percent Cover") +
  scale_color_manual(values = c("mytilus californianus" = "blue", "pisaster ochraceus" = "red")) +
  ggtitle("Mean Percent Cover of Mytilus and Pisaster Over Time") +
  theme_minimal()









