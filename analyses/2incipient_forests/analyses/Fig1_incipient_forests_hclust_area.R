
#Joshua G. Smith; jossmith@mbayaq.org

rm(list=ls())


######
#required packages
librarian::shelf(tidyverse, sf, raster, shiny, tmap)

#set directories 
basedir <- "/Volumes/seaotterdb$/kelp_recovery/data"
localdir <- "/Users/jossmith/Documents/Data/landsat"
figdir <- here::here("analyses","2incipient_forests","figures")
output <- here::here("analyses","2incipient_forests","output")

#read landsat 
final_data <- st_read(file.path(output, "landsat_cluster_area.geojson")) %>%
  mutate(site_numeric = as.numeric(site_num))

#read state
ca_counties <- st_read(file.path(basedir, "gis_data/raw/ca_county_boundaries/s7vc7n.shp")) 

# Get land
usa <- rnaturalearth::ne_states(country="United States of America", returnclass = "sf")
foreign <- rnaturalearth::ne_countries(country=c("Canada", "Mexico"), returnclass = "sf")

################################################################################
#set incipient forests
final_data <- final_data %>%
  #set recovery category
  mutate(Incipient = ifelse(site_name == "Butterfly House"|
                              site_name == "Cannery Row"|
                              site_name == "Carmel River Beach" |
                              site_name == "Coast Guard Pier" |
                              site_name == "Lone Cypress"|
                              site_name =="Mono-Lobo" |
                              site_name == "Pescadero East" |
                              site_name == "Pescadero West" |
                              site_name == "Sea Lion Point"|
                              site_name == "Monastery"|
                              site_name == "Stillwater Cove" | 
                              site_name == "Whaler's Cove","Yes","No"),
         incipient = factor(Incipient, levels = c("Yes","No"))) %>%
  #drop sandy sites
  mutate(site_name = factor(site_name))


################################################################################
#plot clusters

#create grid

landsat_hclust <-final_data %>% filter(year == 2023 & quarter == 3) %>%
  #transform to teale
  st_transform(crs=3310)

r <- terra::rast(landsat_hclust, res=30)  # Build a raster
vr_clust <- terra::rasterize(landsat_hclust, r,"site_name") 


# Create the "Monterey" text label
monterey_label <- data.frame(
  x = c(-121.9, -121.98), # x-coordinate for the upper right corner
  y = c(36.64, 36.54),  # y-coordinate for the upper right corner
  label = c("Monterey \nBay", "Carmel \nBay")
)


# Choose a color palette with distinct colors
custom_palette <- RColorBrewer::brewer.pal(n = 8, name = "Dark2")

# Expand the palette to 28 colors by repeating it
custom_palette <- rep(custom_palette, length.out = 28)

# Set "NA" values to be transparent
custom_palette[is.na(custom_palette)] <- "transparent"

hclust_name <- landsat_hclust %>% group_by(site_name) %>%
  st_transform(crs=4326)%>%
  mutate(centroid = st_centroid(geometry))%>%
  mutate(longitude = st_coordinates(centroid)[, 1], 
         latitude = st_coordinates(centroid)[, 2])


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

p1 <- ggplot() +
  # Add clusters
  tidyterra::geom_spatraster(data = as.factor(vr_clust), na.rm = TRUE) +
  scale_fill_manual(values = custom_palette, na.value = "transparent") +
  geom_sf(data = ca_counties, fill = "gray", color = "gray80") +
  coord_sf(xlim = c(-121.99, -121.88), ylim = c(36.519, 36.645), crs = 4326) +
  labs(title = "", tag = "") +
  # Add landmarks
  geom_text(data = monterey_label, mapping = aes(x = x, y = y, label = label),
            size = 3, fontface = "bold") +
  # Add CA inset
  annotation_custom(grob = g1_inset, 
                    xmin = -122.01, 
                    xmax = -121.96,
                    ymin = 36.625) +
  # Add jittered labels for site_name with boxes
  ggrepel::geom_label_repel(
    data = hclust_name,
    aes(x = longitude, y = latitude, label = site_name),
    box.padding = 0.3,
    point.padding = 0.5,
    force = 18,
    size = 2,
    min.segment.length = 0.1,
    segment.color = "grey50"
  ) +
  #add scale bar
  ggsn::scalebar(x.min = -121.99, x.max = -121.88, 
                 y.min = 36.519, y.max = 36.645,
                 #anchor=c(x=-124.7,y=41),
                 location="bottomright",
                 dist = 2, dist_unit = "km",
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
  theme_bw() +  theme(
    plot.tag.position = c(-0.03, 1),
    axis.title = element_blank()) +
  labs(title = "",
       x="",
       y="")+
  guides(fill = "none") +base_theme

p1

################################################################################
#plot kelp forest trends by cluster

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
final_data <- final_data[order(final_data$site_numeric), ]

p2 <- ggplot(final_data %>% filter(year >= 1990) %>%
             filter(!(site_name=="South lobos")), 
             aes(x = year, y = deviation, group = site_name)) +
  geom_line() +
  geom_point() +
  facet_wrap(~reorder(site_name, site_num), ncol = 4, scales = "fixed") +
  # Heatwave
  annotate(geom="rect", xmin=2013.5, xmax=2016.5, ymin=-Inf, ymax=Inf, fill="red", alpha=0.2) +
  geom_segment(aes(x = -Inf, xend = 2013.5, y = 0, yend = 0), linetype = "solid", color = "black") +
  geom_hline(yintercept = 0, linetype = "dotted", color = "black") +  # Add this line for the dotted line
  labs(x = "Year", y = "Standard deviations from baseline (2000-2013)") +
  scale_x_continuous(limits = c(1990, 2023), breaks = seq(1990, 2023, 5)) +  
  #scale_fill_manual(values = c("navyblue","indianred"))+
  #scale_color_manual(values = c("navyblue","indianred")) +
  theme_bw() + base_theme 
p2



ggpubr::ggarrange(p1,p2, widths = c(0.4,0.6))


ggsave(g0, filename=file.path(figdir, "FigX_landsat_atos_trend_2000-2022.png"), 
       width=7, height=8, units="in", dpi=600)









######################OLD













g1 <- ggplot(final_data, aes(x = year, y = deviation, group = site_name)) +
  geom_line() +
  geom_point() +
  facet_wrap(~reorder(site_name, site_num), ncol = 4, scales = "fixed") +
  # Heatwave
  annotate(geom="rect", xmin=2013.5, xmax=2016.5, ymin=-Inf, ymax=Inf, fill="red", alpha=0.2) +
  geom_segment(aes(x = 1984, xend = 2013.5, y = 0, yend = 0), linetype = "solid", color = "black") +
  geom_hline(yintercept = 0, linetype = "dotted", color = "black") +  # Add this line for the dotted line
  labs(x = "Year", y = "Standard deviations from baseline (2000-2013)") +
  #scale_fill_manual(values = c("navyblue","indianred"))+
  #scale_color_manual(values = c("navyblue","indianred")) +
  theme_bw() + base_theme 
g1


# Export figure
ggsave(g1, filename=file.path(figdir, "FigX_landsat_atos_trend_1984-2022.png"), 
       width=7, height=8, units="in", dpi=600)



g2 <- ggplot(final_data %>% filter(year>2013), aes(x = year, y = deviation, group = site_name, color = Incipient, fill = Incipient)) +
  geom_line() +
  geom_point() +
  facet_wrap(~reorder(site_name, site_num), ncol = 4, scales = "free_y") +
  # Heatwave
  annotate(geom="rect", xmin=2013.5, xmax=2016.5, ymin=-Inf, ymax=Inf, fill="red", alpha=0.2) +
  labs(x = "Year", y = "Standard deviations from baseline (2000-2013)") +
  scale_fill_manual(values = c("navyblue","indianred"))+
  scale_color_manual(values = c("navyblue","indianred"))+
  scale_x_continuous(breaks = seq(2014, 2022, by = 2))+  # Adjust breaks for rounded years
  theme_bw() + base_theme 
g2


# Export figure
ggsave(g2, filename=file.path(figdir, "FigX_landsat_atos_trend_2014-2022.png"), 
       width=7, height=8, units="in", dpi=600)



################################################################################
#misc code


#######attempt to plot forage dives

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

