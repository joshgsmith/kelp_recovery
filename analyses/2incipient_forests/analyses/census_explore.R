
#Joshua G. Smith; jossmith@mbayaq.org

rm(list=ls())


######
#required packages
librarian::shelf(tidyverse, sf, raster, shiny, tmap, terra, tidyterra, ggnewscale)

#set directories 
basedir <- "/Volumes/seaotterdb$/kelp_recovery/data"
figdir <- here::here("analyses","2incipient_forests","figures")
output <- here::here("analyses","2incipient_forests","output")

#read landsat averages with scan area
scan_orig <- st_read(file.path(output, "landsat_scan_area.geojson"))

#read otter scan area
scan_area <- st_read(file.path(basedir, "gis_data/raw/otter_scan_area/TrackingMaps_POLY.shp")) 

#read state
ca_counties_orig <- st_read(file.path(basedir, "gis_data/raw/ca_county_boundaries/s7vc7n.shp")) 

# Get land
usa <- rnaturalearth::ne_states(country="United States of America", returnclass = "sf")
foreign <- rnaturalearth::ne_countries(country=c("Canada", "Mexico"), returnclass = "sf")

# read census data
census_orig <- st_read(file.path("/Volumes/seaotterdb$/RESEARCH FORMS AND FILES/Monterey monthly census/QuarterlyCensus/qc_processed/merged_test_JGS.shp"))

################################################################################
# format scan area

scan_plot_area <- scan_area %>% mutate(Name = factor(Name)) %>%
  rename(site_name = Name)%>%
  mutate(Incipient = ifelse(site_name == "3200"|
                              site_name == "Carmel River Beach"|
                              site_name == "Coast Guard Pier" |
                              site_name == "El Torito" |
                              site_name == "Hopkins Marine Station East"|
                              site_name == "Lone Cypress" |
                              site_name =="Monastery" |
                              site_name == "Monterey Bay Inn" |
                              site_name == "Pescadero Point East" |
                              site_name == "Pescadero Point West"|
                              site_name == "Whaler's Cove","Yes","No"),
         incipient = factor(Incipient, levels = c("Yes","No")))

################################################################################
# process data

#extract date
census_build1 <- census_orig %>%
  mutate(
    year = lubridate::year(date_time),
    month = lubridate::month(date_time),
    day = lubridate::day(date_time),
    #set quarter
    quarter = case_when(
      month %in% c(1, 2, 3) ~ 1,
      month %in% c(4, 5, 6) ~ 2,
      month %in% c(7, 8, 9) ~ 3,
      month %in% c(10, 11, 12) ~ 4
    )
  ) %>% #set geometry
  st_transform(crs = st_crs(scan_plot_area))%>%
  #join with scan area
  st_join(scan_plot_area)

#determine total counts by quarter and year for foraging otters
census_build2 <- census_build1 %>%
                  filter(behavior == "foraging")%>%
                 group_by(year, quarter)%>%
                  dplyr::summarize(n_independ_total = sum(num_indepen)) 

#determine total counts by scan area and quarter for foraging otters
census_build3 <- census_build1 %>%
  filter(behavior == "foraging")%>%
  group_by(year, quarter, site_name)%>%
  dplyr::summarize(n_site_foragers = sum(num_indepen)) %>%
  #join with total
  st_join(census_build2) %>%
  #determine relative foraging
  mutate(rel_forage = n_site_foragers / n_independ_total)

#join with incipient forsts

census_build4 <- scan_plot_area %>% st_join(census_build3)


################################################################################
# plot

census_plot <- census_build4 %>% filter(year.x == 2021)

# Theme
base_theme <-  theme(axis.text=element_text(size=6, color = "black"),
                     axis.title=element_text(size=7,color = "black"),
                     plot.tag=element_text(size=7,color = "black"),
                     plot.title=element_text(size=7,color = "black", face = "bold"),
                     # Gridlines
                     panel.grid.major = element_blank(), 
                     panel.grid.minor = element_blank(),
                     panel.background = element_blank(), 
                     axis.line = element_line(colour = "black"),
                     # Legend
                     legend.key.size = unit(0.3, "cm"), 
                     legend.key = element_rect(fill=alpha('blue', 0)),
                     legend.spacing.y = unit(0.1, "cm"),  
                     legend.text=element_text(size=5,color = "black"),
                     legend.title=element_text(size=6,color = "black"),
                     #legend.key.height = unit(0.1, "cm"),
                     legend.background = element_rect(fill=alpha('blue', 0)),
                     #facets
                     strip.text = element_text(size=6, face = "bold",color = "black", hjust=0),
                     strip.background = element_blank())

p1 <- ggplot() +
  geom_sf(data = census_plot, aes(fill = rel_forage)) +
  ggnewscale::new_scale_fill() +  # Start a new fill scale
  geom_sf(data = scan_plot_area, aes(fill = Incipient), color = "black",
          alpha = 0.4) +
  scale_fill_manual(values = c("lightpink", "transparent"),
                    breaks = c("Yes", "No"),
                    labels = c("", "")) +
  labs(fill = "Incipient \nforest")+
  #add land
  geom_sf(data = ca_counties_orig, fill = "gray", color = "gray80") +
  coord_sf(xlim = c(-121.99, -121.88), ylim = c(36.519, 36.645), crs = 4326)+
  labs(title = "", tag = "")+
  #reorder legend
  guides(
    fill = guide_legend(order = 2),  
    fill_new_scale = guide_legend(order = 3),  
    fill_new_scale1 = guide_legend(order = 1)  
  ) +  
  theme_bw() + base_theme +theme(axis.text.y = element_blank(),
                                 plot.tag.position = c(-0.03, 1),
                                 axis.title=element_blank())
p1











