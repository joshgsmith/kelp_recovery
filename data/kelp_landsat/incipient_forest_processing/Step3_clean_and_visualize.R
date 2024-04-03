
#Joshua G. Smith; jossmith@mbayaq.org

rm(list=ls())


######
#required packages
librarian::shelf(tidyverse, sf, raster, shiny, tmap)

#set directories 
basedir <- "/Volumes/seaotterdb$/kelp_recovery/data"
localdir <- "/Users/jossmith/Documents/Data/landsat"
figdir <- here::here("data","kelp_landsat","incipient_forest_processing","Figures")
output <- here::here("analyses","2incipient_forests","output")

#read landsat 
final_data <- st_read(file.path(output, "area_data3.geojson")) %>%
  mutate(site_numeric = as.numeric(cluster)) %>%
  st_cast("POINT")

#read state
ca_counties <- st_read(file.path(basedir, "gis_data/raw/ca_county_boundaries/s7vc7n.shp")) 

# Get land
usa <- rnaturalearth::ne_states(country="United States of America", returnclass = "sf")
foreign <- rnaturalearth::ne_countries(country=c("Canada", "Mexico"), returnclass = "sf")



###############################################################################
#clean up spatial extent 

xlims <- c(-121.937462, -121.935460, -121.996530, -121.996530)
ylims <- c(36.640507, 36.579086, 36.581153, 36.641310)

box_coords <- tibble(x = xlims, y = ylims) %>% 
  st_as_sf(coords = c("x", "y")) %>% 
  st_set_crs(st_crs(final_data))

#get the bounding box of the two x & y coordintates, make sfc
bounding_box <- st_bbox(box_coords) %>% st_as_sfc()
plot(bounding_box)

# Filter out the points falling within the bounding box
landsat_build1 <- st_difference(final_data, bounding_box) 


#clean up
xlims <- c(-121.891489, -122.010802, -121.991380, -121.891489)
ylims <- c(36.651895, 36.658547, 36.519945, 36.519945)


box_coords <- tibble(x = xlims, y = ylims) %>% 
  st_as_sf(coords = c("x", "y")) %>% 
  st_set_crs(st_crs(landsat_build1))

#get the bounding box of the two x & y coordintates, make sfc
bounding_box <- st_bbox(box_coords) %>% st_as_sfc()

# Filter out the points falling within the bounding box
landsat_build2 <- st_intersection(landsat_build1, bounding_box) 


plot(landsat_build2 %>% filter(year == 2023))

################################################################################
#rename clusters

#rename sites and create table
landsat_build3 <- landsat_build2 %>%
  mutate(
  #set cluster order
  site_num = case_when(
    cluster == 79 ~ 1,
    cluster == 78 ~ 2,
    cluster == 77 ~ 3,
    cluster == 76 ~ 4,
    cluster == 75 ~ 5,
    cluster == 74 ~ 6,
    cluster == 73 ~ 7,
    cluster == 72 ~ 8, 
    cluster == 71 ~ 9,
    cluster == 70 ~ 10,
    cluster == 69 ~ 11,
    cluster == 68 ~ 12,
    cluster == 67 ~ 13,
    cluster == 66 ~ 14,
    cluster == 64 ~ 15,
    cluster == 101 ~ 16,
    cluster == 100 ~ 17,
    cluster == 45 ~ 18,
    cluster == 44 ~ 19,
    cluster == 1 ~ 20,
    cluster == 2 ~ 21,
    cluster == 4 ~ 22,
    cluster == 85 ~ 23,
    cluster == 86 ~ 24,
    cluster == 6 ~ 25,
    cluster == 7 ~ 26,
    cluster == 103 ~ 27,
    cluster == 14 ~ 28,
    cluster == 11 ~ 29,
    cluster == 16 ~ 30,
    cluster == 17 ~ 31,
    cluster ==22 ~ 32,
    cluster == 21 ~ 33,
    cluster == 9 ~ 34,
    cluster == 88 ~ 35,
    cluster == 23 ~ 36,
    cluster == 24 ~ 37,
    cluster == 119 ~ 38,
    cluster == 28 ~ 39,
    cluster == 96 ~ 40,
    cluster == 97 ~ 41,
    cluster == 31 ~ 42,
    cluster == 33 ~ 43,
    cluster == 37 ~ 44,
    cluster == 39 ~ 45,
    cluster == 40 ~ 46,
    cluster == 43 ~ 47,
    cluster == 51 ~ 48,
    cluster == 56 ~ 49,
    cluster == 47 ~ 50,
    cluster == 46 ~ 51,
    cluster == 55 ~ 52,
    cluster == 59 ~ 53,
    cluster == 99 ~ 54,
    cluster == 53 ~ 55,
    cluster == 60 ~ 56,
    cluster == 50 ~ 57,
    cluster == 54 ~ 58,
    cluster == 62 ~ 59,
    cluster == 49 ~ 60,
    cluster == 52 ~ 61,
    cluster == 58 ~ 62,
    cluster == 61 ~ 63,
    cluster == 117 ~ 64,
    cluster == 65 ~ 65,
    cluster == 63 ~ 66,
    cluster == 57 ~ 67,
    cluster == 48 ~ 68,
    cluster == 41 ~ 69,
    cluster == 38 ~ 70,
    cluster == 98 ~ 71,
    cluster == 123 ~ 72,
    cluster == 110 ~ 73,
    TRUE ~ NA)
  ) %>%
  #filter to data extent
  filter(!(is.na(cluster)))




################################################################################
#plot timeseries for each cluster to inspect

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


# Plot 
ggplot(landsat_build3 %>% filter(year > 2013), aes(x = year, y = perc_of_max_3)) +
  geom_point() + 
  geom_smooth(se = TRUE) +
  facet_wrap(~ site_num, scales = "free_y") + theme_bw() + base_theme


################################################################################
#determine center coords (THIS WILL BE IMPORTANT TO SAVE LATER)

cluster_coord <- landsat_build3 %>%
  mutate(longitude = sf::st_coordinates(.)[,1],
         latitude = sf::st_coordinates(.)[,2])%>%
  st_drop_geometry()%>%
  group_by(site_num)%>%
  dplyr::summarize(lat = mean(latitude),
                   long = mean(longitude))%>%
  ungroup() %>%
  st_as_sf(coords = c("long", "lat"), crs = 4326)


#inspect

ggplot(data = cluster_coord) +
  #geom_sf() +  # Plot the spatial data
  geom_sf(data = ca_counties, fill = "gray", color = "gray80") +
  geom_sf_text(aes(label = site_num), size = 4, color = "black")+
  coord_sf(xlim = c(-121.99, -121.88), ylim = c(36.519, 36.645), crs = 4326) 



################################################################################
#label incipient clusters

#rename sites and create table
landsat_build4 <- landsat_build3 %>%
  mutate(
    #set cluster order
    incipient = case_when(
      site_num == 1 ~ "Incipient",
      site_num == 2 ~ "Incipient",
      site_num == 3 ~ "Forest",
      site_num == 4 ~ "Forest",
      site_num == 5 ~ "Barren",
      site_num == 6 ~ "Barren",
      site_num == 7 ~ "Barren",
      site_num == 8 ~ "Barren",
      site_num == 9 ~ "Barren",
      site_num == 10 ~ "Barren",
      site_num == 11 ~ "Barren",
      site_num == 12 ~ "Barren",
      site_num == 13 ~ "Barren",
      site_num == 14 ~ "Barren",
      site_num == 15 ~ "Forest",
      site_num == 16 ~ "Incipient",
      site_num == 17 ~ "Barren",
      site_num == 18 ~ "Barren",
      site_num == 19 ~ "Incipient",
      site_num == 20 ~ "Barren",
      site_num == 21 ~ "Barren",
      site_num == 22 ~ "Barren",
      site_num == 23 ~ "Forest",
      site_num == 24 ~ "Barren",
      site_num == 25 ~ "Forest",
      site_num == 26 ~ "Incipient",
      site_num == 27 ~ "Incipient",
      site_num == 28 ~ "Barren",
      site_num == 29 ~ "Incipient",
      site_num == 30 ~ "Barren",
      site_num == 31 ~ "Barren",
      site_num == 32 ~ "Incipient",
      site_num == 33 ~ "Incipient",
      site_num == 34 ~ "Incipient",
      site_num == 35 ~ "Incipient",
      site_num == 36 ~ "Barren",
      site_num == 37 ~ "Barren",
      site_num == 38 ~ "Barren",
      site_num == 39 ~ "Incipient",
      site_num == 40 ~ "Incipient",
      site_num == 41 ~ "Forest",
      site_num == 42 ~ "Barren",
      site_num == 43 ~ "Incipient",
      site_num == 44 ~ "Incipient",
      site_num == 45 ~ "Incipient",
      site_num == 46 ~ "Incipient",
      site_num == 47 ~ "Barren",
      site_num == 48 ~ "Incipient",
      site_num == 49 ~ "Forest",
      site_num == 50 ~ "Forest",
      site_num == 51 ~ "Barren",
      site_num == 52 ~ "Forest",
      site_num == 53 ~ "Forest",
      site_num == 54 ~ "Barren",
      site_num == 55 ~ "Forest",
      site_num == 56 ~ "Forest",
      site_num == 57 ~ "Barren",
      site_num == 58 ~ "Barren",
      site_num == 59 ~ "Forest",
      site_num == 60 ~ "Forest",
      site_num == 61 ~ "Barren",
      site_num == 62 ~ "Incipient",
      site_num == 63 ~ "Forest",
      site_num == 64 ~ "Forest",
      site_num == 65 ~ "Incipient",
      site_num == 66 ~ "Incipient",
      site_num == 67 ~ "Incipient",
      site_num == 68 ~ "Incipient",
      site_num == 69 ~ "Barren",
      site_num == 70 ~ "Incipient",
      site_num == 71 ~ "Forest",
      site_num == 72 ~ "Barren",
      site_num == 73 ~ "Forest",
      TRUE ~ NA)
  ) %>% filter(year == 2023)


# Plot
ggplot() +
  geom_sf(data = landsat_build4, aes(color = incipient)) +
  scale_color_manual(values = c("Forest" = "forestgreen", "Barren" = "purple", "Incipient" = "orange"), name = "Incipient") +
  coord_sf(crs = 4326) +
  theme_minimal() 



################################################################################
#plot everything

cluster_coord <- landsat_build3 %>%
  mutate(longitude = sf::st_coordinates(.)[,1],
         latitude = sf::st_coordinates(.)[,2])%>%
  st_drop_geometry()%>%
  group_by(site_num)%>%
  dplyr::summarize(lat = mean(latitude),
                   long = mean(longitude))


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
#g1_inset

p1 <- ggplot() +
  # Add clusters
  geom_sf(data = landsat_build4, aes(color = incipient)) +
  scale_color_manual(values = c("Forest" = "forestgreen", "Barren" = "purple", "Incipient" = "orange"), name = "Site type") +
  geom_sf(data = ca_counties, fill = "gray", color = "gray80") +
  labs(title = "", tag = "") +
  # Add landmarks
  #geom_text(data = monterey_label, mapping = aes(x = x, y = y, label = label),
  #        size = 3, fontface = "bold") +
  # Add CA inset
  annotation_custom(grob = g1_inset, 
                    xmin = -122.01, 
                    xmax = -121.96,
                    ymin = 36.625) +
  # Add jittered labels for site_name with boxes
  ggrepel::geom_label_repel(
    data = cluster_coord,
    aes(x = long, y = lat, label = site_num),
    box.padding = 0.3,
    point.padding = 0.5,
    force = 18,
    size = 2,
    min.segment.length = 0.1,
    segment.color = "black"
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
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank())+
  #guides(fill = guide_legend(override.aes = list(size = 3))) +
  base_theme+
  coord_sf(xlim = c(-121.99, -121.88), ylim = c(36.519, 36.645), crs = 4326) 

p1


ggsave(p1,  filename=file.path(figdir, "Cluster_site_type_map.png"), width = 7.5, height = 9.5, units = "in",
       bg = "white", dpi = 600)




# Plot 
###remove year filter on landsat_buld4 before running
p2 <- ggplot(landsat_build4 %>% filter(year > 2013), aes(x = year, y = perc_of_max_3,
                                                         color = incipient)) +
  geom_point() + 
  geom_smooth(se = TRUE) +
  scale_color_manual(values = c("Forest" = "forestgreen", "Barren" = "purple", "Incipient" = "orange"), name = "Site type") +
  facet_wrap(~ site_num, scales = "free_y") + theme_bw() + base_theme
p2

ggsave(p2,  filename=file.path(figdir, "Cluster_timeseries.png"), width = 16, height = 10, units = "in",
       bg = "white", dpi = 600)



