
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
landsat_dat <- st_read(file.path(basedir, "kelp_landsat/processed/monterey_peninsula/landsat_mpen_1984_2023_points_withNAs.shp"))%>%
  #filter data to reduce memory load
  filter(year > 1999) %>% filter(quarter == 3)

# Get land
ca_counties_orig <- st_read(file.path(basedir, "gis_data/raw/ca_county_boundaries/s7vc7n.shp")) 

################################################################################
#
# Step 1: Calculate baseline average for each point
baseline_avg <- landsat_dat %>%
  filter(year >= 2000 & year <= 2013) %>%
  group_by(geometry) %>%
  summarize(baseline_avg = mean(biomass),
            baseline_sd = sd(biomass))

# Step 2: Calculate standard deviation from observed value in 2020 versus baseline average
result <- landsat_dat %>%
  st_join(baseline_avg) %>%
  mutate(std_dev = (biomass - baseline_avg) / baseline_sd) %>%
  #clean up
  mutate(std_dev = ifelse(std_dev %in% c("NaN","Inf"),NA,std_dev)) %>%
  filter(year == 2022 & quarter == 3) %>%
  #get only positive biomass
  filter(biomass > 0)



#Step raster for focal year and quarter

# transform landsat data to Teale Albers
rast_build1 <- st_transform(result, crs = 3310) %>% 
  mutate(biomass = ifelse(biomass == 0, NA, biomass))  %>% filter (year == 2022 & quarter ==3)  %>% filter(std_dev <2)

#define blank raster
r <- rast(rast_build1, res=30)

landsat_rast_2022 <- rasterize(rast_build1, r, field = "std_dev") 


################################################################################
#plot

ggplot() +
  #add kelp extent
  tidyterra::geom_spatraster(data = landsat_rast_2022, na.rm = TRUE) +
  tidyterra::scale_fill_whitebox_c(
    palette = "atlas",
    na.value = NA
  ) +
  labs(fill = expression("Kelp biomass" ~(kg~"/"~"30m"^{-2})))+
  # Plot land
  #add land
  geom_sf(data = ca_counties_orig, fill = "gray", color = "gray80") +
  coord_sf(xlim = c(-121.99, -121.88), ylim = c(36.519, 36.645), crs = 4326)


