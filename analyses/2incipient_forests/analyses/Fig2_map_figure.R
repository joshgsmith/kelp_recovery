

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

#read landsat raw
landsat_orig <- st_read(file.path(basedir, "kelp_landsat/processed/monterey_peninsula/landsat_mpen_1984_2023_points_withNAs.shp"))

#read foraging data
forage_orig <- read.csv(file.path(basedir, "foraging_data/processed/foraging_data_2016_2023.csv")) 

# Get land
usa <- rnaturalearth::ne_states(country="United States of America", returnclass = "sf")
foreign <- rnaturalearth::ne_countries(country=c("Canada", "Mexico"), returnclass = "sf")


################################################################################
#Step 1 - raster for focal year and quarter

# transform landsat data to Teale Albers
rast_build1 <- st_transform(landsat_orig, crs = 3310) %>% 
  mutate(biomass = ifelse(biomass == 0, NA, biomass))  %>% filter (year == 2022 & quarter ==3) 

#define blank raster
r <- rast(rast_build1, res=30)

landsat_rast_2022 <- rasterize(rast_build1, r, field = "biomass", fun = mean)



rast_build2 <- st_transform(landsat_orig, crs = 3310) %>% 
  mutate(biomass = ifelse(biomass == 0, NA, biomass))  %>% filter (year == 2017 & quarter ==3) 
landsat_rast_2017 <- rasterize(rast_build2, r, field = "biomass", fun = mean)

################################################################################
#Step 2 - define max kelp extent

#select MPEN as focal region 
plot_dat <- rast_build1 %>% filter (latitude >= 36.510140 &
                                       latitude <= 36.670574) %>% st_transform(crs=3310)

#define 0 cover as historical kelp footprint
na_dat <- landsat_orig %>% filter (latitude >= 36.510140 &
                                    latitude <= 36.670574, 
                                  biomass == 0)

#transform landsat data to Teale Albers
t_dat <- st_transform(plot_dat, crs=3310) %>% 
  mutate(biomass = ifelse(biomass==0,NA,biomass))

kelp_historic <- st_transform(na_dat, crs=3310) # %>% mutate(biomass = ifelse(biomass==0,1,NA))

#create grid
r <- rast(t_dat, res=30)  # Builds a blank raster of given dimensions and resolution  
vr <- rasterize(t_dat, r,"biomass", resolution = 30) 

kelp_na <- rasterize(kelp_historic, r,"biomass", resolution = 30) 



################################################################################
#Step 3 - process scan area

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
#Step 4 - plot

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

# Build inset
g1_inset <-  ggplotGrob(
  ggplot() +
    # Plot land
    geom_sf(data=foreign, fill="grey80", color="white", lwd=0.3) +
    geom_sf(data=usa, fill="grey80", color="white", lwd=0.3) +
    # Plot box
    annotate("rect", xmin=-122.6, xmax=-121, ymin=36.2, ymax=37.1, color="black", fill=NA, lwd=0.6) +
    # Label regions
    #geom_text(data=region_labels, mapping=aes(y=lat_dd, label=region), x= -124.4, hjust=0, size=2) +
    # Labels
    labs(x="", y="") +
    # Crop
    coord_sf(xlim = c(-124.5, -117), ylim = c(32.5, 42)) +
    # Theme
    theme_bw() + base_theme +
    theme( plot.margin = unit(rep(0, 4), "null"),
           panel.margin = unit(rep(0, 4), "null"),
           panel.background = element_rect(fill='transparent'), #transparent panel bg
           # plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
           axis.ticks = element_blank(),
           axis.ticks.length = unit(0, "null"),
           axis.ticks.margin = unit(0, "null"),
           axis.text = element_blank(),
           axis.title=element_blank(),
           axis.text.y = element_blank())
)
g1_inset

# Create the "Monterey" text label
monterey_label <- data.frame(
  x = c(-121.9, -121.98), # x-coordinate for the upper right corner
  y = c(36.64, 36.54),  # y-coordinate for the upper right corner
  label = c("Monterey \nBay", "Carmel \nBay")
)

 
p1 <- ggplot() +
  geom_sf(data = scan_plot_area, aes(fill = Incipient), color = NA,
          alpha = 0.4) +
  scale_fill_manual(values = c("lightpink", "transparent"),
                    breaks = c("Yes", "No"),
                    labels = c("", "")) +
  labs(fill = "Incipient \nforest")+
  ggnewscale::new_scale_fill()+
  #add historic kelp extent
  tidyterra::geom_spatraster(data = kelp_na, na.rm = TRUE) +
  tidyterra::scale_fill_hypso_c(
    palette = "moon_bathy",
    na.value = NA,
    labels = NULL
  ) +
  guides(fill = guide_legend(override.aes = list(size = 3))) +
  labs(fill = "Max kelp \nextent")+
  #add observed landsat for 2022 Q3
  ggnewscale::new_scale_fill()+
  tidyterra::geom_spatraster(data = landsat_rast_2022, na.rm = TRUE) +
  tidyterra::scale_fill_whitebox_c(
    palette = "atlas",
    na.value = NA
  ) +
  labs(fill = expression("Kelp biomass" ~(kg~"/"~"30m"^{-2})))+
  #labs(fill = expression(atop("Kelp biomass", ~(kg~"per"~"30m"^{-2}))))+
  #add scan area
  ggnewscale::new_scale_fill() +  # Start a new fill scale
  geom_sf(data = scan_plot_area, fill = NA, color = "gray20") +
  #add land
  geom_sf(data = ca_counties_orig, fill = "gray", color = "gray80") +
  coord_sf(xlim = c(-121.99, -121.88), ylim = c(36.519, 36.645), crs = 4326)+
  labs(title = "2022", tag = "B")+
  #reorder legend
  guides(
    fill = guide_legend(order = 2),  
    fill_new_scale = guide_legend(order = 3),  
    fill_new_scale1 = guide_legend(order = 1)  
  ) +  
  #add landmarks
  geom_text(data=monterey_label, mapping=aes(x=x, y=y, label=label),
            size=3, fontface= "bold") +
  theme_bw() + base_theme +theme(axis.text.y = element_blank(),
                                      plot.tag.position = c(-0.03, 1),
                                 axis.title=element_blank())
p1


p2 <- ggplot() +
  #add historic kelp extent
  tidyterra::geom_spatraster(data = kelp_na, na.rm = TRUE) +
  tidyterra::scale_fill_hypso_c(
    palette = "moon_bathy",
    na.value = NA,
    labels = NULL
  ) +
  guides(fill = guide_legend(override.aes = list(size = 3))) +
  labs(fill = "Max kelp \nextent")+
  #add observed landsat for 2022 Q3
  ggnewscale::new_scale_fill()+
  tidyterra::geom_spatraster(data = landsat_rast_2017, na.rm = TRUE) +
  tidyterra::scale_fill_whitebox_c(
    palette = "atlas",
    na.value = NA
  ) +
  labs(fill = expression(atop("Kelp biomass", ~(kg~"per"~"30m"^{-2}))))+
  #add scan area
  ggnewscale::new_scale_fill() +  # Start a new fill scale
  geom_sf(data = scan_plot_area, fill = NA, color = "gray20") +
  #add land
  geom_sf(data = ca_counties_orig, fill = "gray", color = "gray80") +
  #add inset
  # Add CA inset
  annotation_custom(grob = g1_inset, 
                    xmin = -122.01, 
                    xmax = -121.96,
                    ymin = 36.625) +
  #add scale bar
  ggsn::scalebar(x.min = -121.99, x.max = -121.88, 
                y.min = 36.519, y.max = 36.645,
               #anchor=c(x=-124.7,y=41),
              location="bottomright",
             dist = 2.5, dist_unit = "km",
            transform=TRUE, 
           model = "WGS84",
          st.dist=0.02,
         st.size=2,
        border.size=.5,
     height=.02
 )+
  #add north arrow
  ggsn::north(x.min = -121.99, x.max = -121.88, 
              y.min = 36.519, y.max = 36.65,
              location = "topright", 
              scale = 0.05, 
              symbol = 10)+
  #add landmarks
  geom_text(data=monterey_label, mapping=aes(x=x, y=y, label=label),
            size=3, fontface= "bold") +
  coord_sf(xlim = c(-121.99, -121.88), ylim = c(36.519, 36.645), crs = 4326)+
  labs(title = "2017", tag = "A")+
  theme_bw() + base_theme + theme(legend.position = "none", plot.tag.position = c(0.05, 1), axis.title=element_blank())
p2


#make tile plot

tile_theme <-  theme(axis.text=element_text(size=6, color = "black"),
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
                     strip.text = element_text(size=6, face = "bold",color = "black"),
                     strip.background = element_blank())


# Convert 'year' to numeric
scan_orig$year <- as.factor(as.character(scan_orig$year))

#set Incipient order

scan_orig$Incipient <- factor(scan_orig$Incipient, levels = c("Yes","No")) 


# Create a tile plot
p3 <- ggplot(scan_orig %>% mutate(Incipient = ifelse(
                                Incipient == "Yes", "Incipient","Nonincipient")), aes(x = year, y = reorder(site_order,-site_order), fill = deviation)) +
  geom_tile() +
  facet_grid(Incipient~., scales = "free_y",space="free_y") +
  #scale_fill_gradient2(low = "navyblue",mid="gray80",high="gray80")+ 
  #scale_fill_gradientn(colors = custom_colors, values = scales::rescale(c(0, 1))) +
  scale_fill_gradientn(colors = c("forestgreen", "#FFCC80", "#FFCC80"),
                        values = scales::rescale(c(-2.4, 0, 0.1, 3.1))) +
  labs(x = "Year", y = "Site", fill = "Deviation", tag = "C") +
  #plot marine heatwave
  geom_vline(xintercept = which(levels(factor(scan_orig$year)) %in% c("2014", "2016")), 
             linetype = "dashed", color = "black") +
  #geom_vline(xintercept = 2017.5, linetype = "solid", color = "black")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + # Rotate x-axis labels for better readability
  theme_bw() + tile_theme + theme(axis.text.y = element_blank())+
  scale_x_discrete(breaks=seq(2000, 2022, 3)) 
p3


p <- gridExtra::grid.arrange(p2, p1, nrow=1, widths = c(0.436,0.574))
p

# Export figure
#ggsave(p, filename=file.path(figdir, "Fig2_map_figure.png"), 
 #      width=7, height=5, units="in", dpi=600)


# Adjust the width of p3
p3_adjusted <- gridExtra::arrangeGrob(p3, widths = unit(0.6, "npc"))

# Arrange p and the adjusted p3 in two rows
p_final <- gridExtra::arrangeGrob(p, p3_adjusted, ncol = 1, heights = c(0.6, 0.4))

# Save the final figure
ggsave(p_final, filename = file.path(figdir, "Fig2_map_figure_new.png"), 
       width = 7, height = 7, units = "in", dpi = 600)








###plot time series


#####test


# Convert 'year' to numeric
scan_orig$year <- as.numeric(as.character(scan_orig$year))

# Calculate the mean deviation for each year and Incipient group
mean_deviations_post <- scan_orig %>%
  filter(year >= 2018) %>%
  group_by(year, Incipient) %>%
  summarize(mean_deviation = median(deviation, na.rm = TRUE)) %>%
  st_drop_geometry()

mean_deviations_pre <- scan_orig %>%
  filter(year < 2018) %>%
  group_by(year) %>%
  summarize(mean_deviation = median(deviation, na.rm = TRUE))%>%
  st_drop_geometry()

# Get the last year in mean_deviations_pre
last_year_pre <- max(mean_deviations_pre$year)

# Get the 2017 deviation value
deviation_2017 <- mean_deviations_pre %>% filter(year == 2017) %>% pull(mean_deviation)

# Add rows to mean_deviations_post for 2017 with both Incipient values
mean_deviations_post <- rbind(
  mean_deviations_post,
  data.frame(year = 2017, Incipient = "Yes", mean_deviation = deviation_2017),
  data.frame(year = 2017, Incipient = "No", mean_deviation = deviation_2017)
)

# Create a plot
plot <- ggplot() +
  geom_point(data = scan_orig, aes(x = year, y = deviation), color = "black") +
  #plot pre 2017
  geom_line(data = mean_deviations_pre, aes(x = year, y = mean_deviation), size = 1) +
  # Add a point for the last year in mean_deviations_pre
  geom_point(data = mean_deviations_pre %>% filter(year == last_year_pre), aes(x = year, y = mean_deviation), color = "black") +
  #plot 2017 and beyond
  geom_line(data = mean_deviations_post, aes(x = year, y = mean_deviation, color = Incipient, group = Incipient), size = 1) +
  scale_color_manual(values = c("Yes" = "red", "No" = "blue")) +
  labs(x = "Year", y = "Deviation") +
  theme_minimal()

print(plot)


###make tile plot





