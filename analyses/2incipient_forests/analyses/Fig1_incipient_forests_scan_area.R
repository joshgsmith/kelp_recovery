
#Joshua G. Smith; jossmith@mbayaq.org

rm(list=ls())


######
#required packages
librarian::shelf(tidyverse, sf, raster, shiny, tmap)

#set directories 
basedir <- "/Volumes/seaotterdb$/kelp_recovery/data"
figdir <- here::here("analyses","2incipient_forests","figures")
output <- here::here("analyses","2incipient_forests","output")

#read landsat 
final_data <- st_read(file.path(output, "landsat_scan_area.geojson"))

#read state
ca_counties <- st_read(file.path(basedir, "gis_data/raw/ca_county_boundaries/s7vc7n.shp")) 



################################################################################
#plot kelp forest trends by scan area

# Theme
base_theme <-  theme(axis.text=element_text(size=7, color = "black"),
                     axis.title=element_text(size=8,color = "black"),
                     legend.text=element_text(size=7,color = "black"),
                     legend.title=element_text(size=8,color = "black"),
                     plot.tag=element_text(size=8,color = "black"),
                     # Gridlines
                     panel.grid.major = element_blank(), 
                     panel.grid.minor = element_blank(),
                     panel.background = element_blank(), 
                     axis.line = element_line(colour = "black"),
                     # Legend
                     legend.key = element_rect(fill=alpha('blue', 0)),
                     legend.background = element_rect(fill=alpha('blue', 0)),
                     #facets
                     strip.text = element_text(size=6, face = "bold",color = "black", hjust=0),
                     strip.background = element_blank())

# Sort the dataframe by site_order
final_data <- final_data[order(final_data$site_order), ]

g1 <- ggplot(final_data, aes(x = year, y = deviation, group = site_name, color = Incipient, fill = Incipient)) +
  geom_line() +
  geom_point() +
  facet_wrap(~reorder(site_name, site_order), ncol = 4, scales = "fixed") +
  # Heatwave
  annotate(geom="rect", xmin=2013.5, xmax=2016.5, ymin=-Inf, ymax=Inf, fill="red", alpha=0.2) +
  geom_segment(aes(x = 2000, xend = 2013.5, y = 0, yend = 0), linetype = "solid", color = "black") +
  geom_hline(yintercept = 0, linetype = "dotted", color = "black") +  # Add this line for the dotted line
  labs(x = "Year", y = "Standard deviations from baseline (2000-2013)") +
  theme_bw() + base_theme +
  scale_fill_manual(values = c("navyblue","indianred"))+
  scale_color_manual(values = c("navyblue","indianred"))
g1


# Export figure
ggsave(g1, filename=file.path(figdir, "FigX_landsat_sd_trend_2000-2022.png"), 
       width=7, height=8, units="in", dpi=600)





g2 <- ggplot(final_data %>% filter(year>2013), aes(x = year, y = deviation, group = site_name, color = Incipient, fill = Incipient)) +
  geom_line() +
  geom_point() +
  facet_wrap(~reorder(site_name, site_order), ncol = 4, scales = "free_y") +
  # Heatwave
  annotate(geom="rect", xmin=2013.5, xmax=2016.5, ymin=-Inf, ymax=Inf, fill="red", alpha=0.2) +
  labs(x = "Year", y = "Standard deviations from baseline (2000-2013)", title = "Kelp canopy deviations from 2000-2013 baseline") +
  theme_bw() + base_theme +
  scale_fill_manual(values = c("navyblue","indianred"))+
  scale_color_manual(values = c("navyblue","indianred"))+
  scale_x_continuous(breaks = seq(2014, 2022, by = 2))  # Adjust breaks for rounded years
g2

# Export figure
ggsave(g2, filename=file.path(figdir, "FigX_landsat_sd_trend_2014-2022.png"), 
       width=7, height=10, units="in", dpi=600)





###attempt to plot forage dives


# Drop geometry and select relevant columns from forage_build6
forage_build6_nonspatial <- forage_build5 %>%
  st_drop_geometry() %>%
  dplyr::select(year = year.x, site_name, prop_dives)

# Merge final_data and forage_build6_nonspatial based on year and site_name
joined_data <- final_data %>%
  filter(year > 2013) %>%
  left_join(forage_build6_nonspatial, by = c("year", "site_name"))



# Create the plot
g3 <- ggplot(joined_data, aes(x = year)) +
  geom_line(aes(y = deviation, color = Incipient)) +
  geom_point(aes(y = deviation, color = Incipient, size = coalesce(prop_dives, 0))) +
  facet_wrap(~reorder(site_name, site_order), ncol = 4, scales = "free_y") +
  # Heatwave
  annotate(geom = "rect", xmin = 2013.5, xmax = 2016.5, ymin = -Inf, ymax = Inf, fill = "red", alpha = 0.2) +
  labs(x = "Year") +
  theme_bw() + base_theme +
  scale_color_manual(values = c("navyblue", "indianred")) +
  scale_size_continuous(range = c(1, 6), name = "Relative dive \nfrequency") +
  scale_y_continuous(name = "Kelp biomass (standard deviations from 2000-2013)") +
  scale_x_continuous(breaks = seq(2014, 2022, by = 2))  # Adjust breaks for rounded years

g3


# Export figure
ggsave(g3, filename=file.path(figdir, "FigX_landsat_foraging_trend_2014-2022.png"), 
       width=7.5, height=8, units="in", dpi=600)

