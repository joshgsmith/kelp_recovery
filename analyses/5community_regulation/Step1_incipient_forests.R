
#Joshua G. Smith; jossmith@mbayaq.org

rm(list=ls())


######
#required packages
librarian::shelf(tidyverse, sf, raster, shiny, tmap)

#set directories 
basedir <- "/Volumes/seaotterdb$/kelp_recovery/data"
figdir <- here::here("analyses","5community_regulation","figures")

#read landsat dat
landsat_dat <- st_read(file.path(basedir, "kelp_landsat/processed/monterey_peninsula/landsat_mpen_1984_2023_points_withNAs.shp"))%>%
                  #filter data to reduce memory load
                  filter(year > 1999)

#read state
ca_counties <- st_read(file.path(basedir, "gis_data/raw/ca_county_boundaries/s7vc7n.shp")) 

#read otter scan area
scan_area <- st_read(file.path(basedir, "gis_data/raw/otter_scan_area/TrackingMaps_POLY.shp")) 


################################################################################
#Step 2 - scale biomass to scan area polygon extent 

scan_area <- scan_area %>% rename(site_name = Name) %>% st_transform(crs = st_crs(landsat_dat)) %>%
                #set order upcoast to downcoast
                mutate(site_order = case_when(
                  site_name == "Costco Beach" ~ 1,
                  site_name == "Monterey Beach Hotel (Tides Hotel)" ~ 2,
                  site_name == "Condos" ~ 3,
                  site_name == "Naval Post Graduate School" ~ 4,
                  site_name == "Del Monte" ~5,
                  site_name == "Monterey Harbor" ~6,
                  site_name == "Coast Guard Pier" ~7,
                  site_name == "Monterey Bay Inn" ~8,
                  site_name == "El Torito" ~9,
                  site_name == "MBA" ~10,
                  site_name == "Hopkins Marine Station East" ~11,
                  site_name == "Hopkins Marine Station West" ~12,
                  site_name == "Berwick" ~13,
                  site_name == "Lover's Point East" ~14,
                  site_name == "Lover's Point West" ~15,
                  site_name == "Otter Point East" ~16, 
                  site_name == "Otter Point West" ~17, 
                  site_name == "Esplanade" ~18,
                  site_name == "Lucas Point" ~19,
                  site_name == "Aumentos" ~20,
                  site_name == "Point Pinos East" ~21,
                  site_name == "Point Pinos West" ~22, 
                  site_name == "Great Tidepool" ~23,
                  site_name == "Asilomar" ~24, 
                  site_name == "Spanish Bay" ~25,
                  site_name == "Moss Beach" ~26,
                  site_name == "Point Joe East" ~27,
                  site_name == "Point Joe West" ~28,
                  site_name == "China Rock" ~29,
                  site_name == "Bird Rock North" ~ 30,
                  site_name == "Bird Rock South" ~31,
                  site_name == "Fanshell" ~32,
                  site_name == "Cypress Point Golf Course" ~33,
                  site_name == "Cypress Point Parking" ~34,
                  site_name == "3200"~35,
                  site_name == "Lone Cypress" ~36,
                  site_name == "Pescadero Point West" ~37,
                  site_name == "Pescadero Point East" ~38,
                  site_name == "Stillwater Cove" ~39,
                  site_name == "Arrowhead Point" ~40,
                  site_name == "Carmel Beach North" ~41,
                  site_name == "Scenic Road" ~42,
                  site_name == "Carmel River Beach"~43,
                  site_name == "Monastery" ~44,
                  site_name == "Whaler's Cove" ~45,
                  site_name == "Sea Lion Point" ~46,
                  TRUE ~ NA
                ))

#inspect

# Calculate centroids of polygons
scan_area$centroid <- st_centroid(scan_area)

# Create a data frame for labels and lines
label_data <- data.frame(
  site_name = scan_area$site_name,
  x = st_coordinates(scan_area$centroid)[, "X"],
  y = st_coordinates(scan_area$centroid)[, "Y"],
  xend = st_coordinates(scan_area$centroid)[, "X"],
  yend = st_coordinates(scan_area$centroid)[, "Y"] + 0.002
)

# Scatter plot with jittered site names and lines connecting labels
ggplot(scan_area) +
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
joined_data <- st_join(landsat_dat, scan_area)

# Group by year, quarter, and site_name, and summarize the biomass values
summarized_data <- joined_data %>%
  group_by(year, quarter, site_name, site_order) %>%
  summarize(total_biomass = sum(biomass, na.rm = TRUE))

  
################################################################################
#prep for plot

#determine baseline average kelp biomass
baseline_average <- summarized_data %>%
  #filter(year >= 2007 & year <= 2013) %>%
  filter(year >= 2000 & year <= 2014) %>% #this sets the baseline period. #use year >= 2016 & year <= 2017
  filter(quarter == 3)%>% #filter to quarter 3 for max kelp extent
  group_by(site_name, site_order) %>%
  summarize(baseline_avg = mean(total_biomass))


observed_means <- summarized_data %>%
  filter(quarter == 3)%>%
  #filter(site_name %in% c("Point Pinos East", "Point Pinos West", "Pescadero Point West", "Pescadero Point East")) %>%
  group_by(site_name,site_order, year) %>%
  summarize(observed_avg = mean(total_biomass))


# Assuming baseline_average and observed_means are data frames
final_data <- observed_means %>%
  filter(year >= 2000)%>%
  st_join(baseline_average) %>%
  rename(site_name = site_name.x,
         site_order = site_order.x) %>%
  dplyr::select(-site_name.y, -site_order.y) %>%
  mutate(deviation = (observed_avg - baseline_avg) / sd(observed_avg)) %>%
  filter(!(is.na(site_name))) %>%
  #set recovery category
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
         incipient = factor(Incipient, levels = c("Yes","No"))) %>%
  #drop sandy sites
  mutate(site_name = factor(site_name))%>%
  filter(!(site_name %in% c("Condos", "Naval Post Graduate School", "Del Monte")))

################################################################################
#prep for plot

# Theme
base_theme <-  theme(axis.text=element_text(size=7),
                     axis.title=element_text(size=8),
                     legend.text=element_text(size=7),
                     legend.title=element_text(size=8),
                     plot.tag=element_text(size=8),
                     # Gridlines
                     panel.grid.major = element_blank(), 
                     panel.grid.minor = element_blank(),
                     panel.background = element_blank(), 
                     axis.line = element_line(colour = "black"),
                     # Legend
                     legend.key = element_rect(fill=alpha('blue', 0)),
                     legend.background = element_rect(fill=alpha('blue', 0)),
                     #facets
                     strip.text = element_text(size=7, face = "bold"),
                     strip.background = element_blank())

# Sort the dataframe by site_order
final_data <- final_data[order(final_data$site_order), ]

g1 <- ggplot(final_data, aes(x = year, y = deviation, group = site_name, color = Incipient, fill = Incipient)) +
  geom_line() +
  geom_point() +
  facet_wrap(~reorder(site_name, site_order), ncol = 4, scales = "free_y") +
  # Heatwave
  annotate(geom="rect", xmin=2013.5, xmax=2016.5, ymin=-Inf, ymax=Inf, fill="red", alpha=0.2) +
  geom_segment(aes(x = 2000, xend = 2013.5, y = 0, yend = 0), linetype = "solid", color = "black") +
  geom_hline(yintercept = 0, linetype = "dotted", color = "black") +  # Add this line for the dotted line
  labs(x = "Year", y = "Standard deviations from baseline (2000-2013)", title = "Kelp canopy deviations from 2000-2013 baseline") +
  theme_bw() + base_theme +
  scale_fill_manual(values = c("navyblue","indianred"))+
  scale_color_manual(values = c("navyblue","indianred"))
g1



# Export figure
ggsave(g1, filename=file.path(figdir, "FigX_landsat_sd_trend_2000-2022.png"), 
       width=7, height=10, units="in", dpi=600)


