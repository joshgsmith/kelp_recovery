
#Joshua G. Smith; jossmith@mbayaq.org

rm(list=ls())


######
#required packages
librarian::shelf(tidyverse, sf, raster, shiny, tmap)

#set directories 
basedir <- "/Volumes/seaotterdb$/kelp_recovery/data"
figdir <- here::here("analyses","2forage_behavior","figures")

#read landsat dat
landsat_dat <- st_read(file.path(basedir, "kelp_landsat/processed/monterey_peninsula/landsat_mpen_1984_2023_points_withNAs.shp"))

#read state
ca_counties <- st_read(file.path(basedir, "gis_data/raw/ca_county_boundaries/s7vc7n.shp")) 

#read otter scan area
scan_area <- st_read(file.path(basedir, "gis_data/raw/otter_scan_area/TrackingMaps_POLY.shp")) 


################################################################################
#Step 2 - scale biomass to scan area polygon extent 

scan_area <- scan_area %>% rename(site_name = Name) 

# Join the point data with the polygon data based on spatial location
joined_data <- st_join(rast_build1, scan_area)

# Group by year, quarter, and site_name, and summarize the biomass values
summarized_data <- joined_data %>%
  group_by(year, quarter, site_name) %>%
  summarize(total_biomass = sum(biomass, na.rm = TRUE))

# Print the summarized data
print(summarized_data)

  
################################################################################
#prep for plot

baseline_average <- summarized_data %>%
  #filter(year >= 2007 & year <= 2013) %>%
  filter(year >= 2016 & year <= 2017) %>%
  filter(quarter == 3)%>%
  group_by(site_name) %>%
  summarize(baseline_avg = mean(total_biomass))

observed_means <- summarized_data %>%
  filter(quarter == 3)%>%
  #filter(site_name %in% c("Point Pinos East", "Point Pinos West", "Pescadero Point West", "Pescadero Point East")) %>%
  group_by(site_name, year) %>%
  summarize(observed_avg = mean(total_biomass))

# Assuming baseline_average and observed_means are data frames
final_data <- observed_means %>%
  filter(year >= 2014)%>%
  st_join(baseline_average, by = "site_name") %>%
  mutate(deviation = (observed_avg - baseline_avg) / sd(observed_avg)) %>%
  filter(!(is.na(site_name.x))) %>%
  #set recovery category
  mutate(Incipient = ifelse(site_name.x == "3200"|
                              site_name.x == "Carmel River Beach"|
                              site_name.x == "Coast Guard Pier" |
                              site_name.x == "El Torito" |
                              site_name.x == "Hopkins Marine Station East"|
                              site_name.x == "Lone Cypress" |
                              site_name.x =="Monastery" |
                              site_name.x == "Monterey Bay Inn" |
                              site_name.x == "Pescadero Point East" |
                              site_name.x == "Pescadero Point West"|
                              site_name.x == "Whaler's Cove","Yes","No"),
         incipient = factor(Incipient, levels = c("Yes","No")))


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



g1 <- ggplot(final_data, aes(x = year, y = deviation, group = site_name.x, color = Incipient, fill = Incipient)) +
  geom_line() +
  geom_point() +
  facet_wrap(~site_name.x, ncol = 5, scales = "fixed") +
  geom_hline(yintercept = 0, linetype = "dotted", color = "black") +  # Add this line for the dotted line
  labs(x = "Year", y = "Deviation (Standard Deviations from 2016-17)", title = "Observed Annual Mean vs. Baseline Deviation by Site") +
  theme_bw() + base_theme +
  scale_fill_manual(values = c("navyblue","indianred"))+
  scale_color_manual(values = c("navyblue","indianred"))
g1



# Export figure
ggsave(g1, filename=file.path(figdir, "FigX_landsat_sd_trend.png"), 
       width=9.5, height=10, units="in", dpi=600)



