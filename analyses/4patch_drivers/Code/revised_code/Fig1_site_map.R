#Communtiy change analyses
#Joshua G. Smith; jossmith@mbayaq.org

rm(list=ls())

librarian::shelf(tidyverse, here, tidyverse, tidync, sf, raster, terra, ggtext, ggsn)



################################################################################
#set directories and load data
basedir <- "/Volumes/seaotterdb$/kelp_recovery/"
figdir <- here::here("analyses","4patch_drivers","Figures")

#load standardized dat
stan_dat <- read.csv(file.path(basedir, "data/subtidal_monitoring/processed/kelp_stan_CC.csv")) 

#read state
ca_counties <- st_read(file.path(basedir, "data/gis_data/raw/ca_county_boundaries/s7vc7n.shp")) 

# Get land
usa <- rnaturalearth::ne_states(country="United States of America", returnclass = "sf")
foreign <- rnaturalearth::ne_countries(country=c("Canada", "Mexico"), returnclass = "sf")


################################################################################
#Identify distinct sites

site_locations_sf <- stan_dat %>% dplyr::select(site, latitude, longitude) %>% distinct() %>%
                    #make sf
                    st_as_sf(coords = c("longitude", "latitude"), crs = 4326)

site_locations_labels <- stan_dat %>% dplyr::select(site, latitude, longitude) %>% distinct() %>%
                                mutate(site = str_to_title(gsub("_", " ", site)),
                                       site = str_replace(site, "Dc", "DC"),
                                       site = str_replace(site, "Uc", "UC"))
                                       


################################################################################
#step 1 - build base map

# create a bounding box with the desired latitude and longitude values
bbox <- st_bbox(c(xmin = -122.08, ymin = 36.4, xmax = -121.7, ymax = 36.670574), crs = st_crs(ca_counties))

# filter the ca_counties object using the bounding box
ca_counties_mpen <- st_intersection(ca_counties, st_as_sfc(bbox))

# landmakrs
landmarks <- tibble(place=c("Monterey \nPeninsula"),
                 long=c(-121.927984),
                 lat=c(36.585))



my_theme <-  theme(axis.text=element_text(size=6),
                   axis.text.y = element_text(angle = 90, hjust = 0.5),
                   axis.title=element_text(size=8),
                   plot.tag=element_text(size=8, face="bold"),
                   plot.title =element_text(size=7, face="bold"),
                   # Gridlines 
                   panel.grid.major = element_blank(), 
                   panel.grid.minor = element_blank(),
                   #panel.background = element_blank(), 
                   axis.line = element_line(colour = "black"),
                   # Legend
                   legend.key = element_blank(),
                   legend.background = element_rect(fill=alpha('blue', 0)),
                   legend.key.height = unit(1, "lines"), 
                   legend.text = element_text(size = 6),
                   legend.title = element_text(size = 7),
                   #legend.spacing.y = unit(0.75, "cm"),
                   #facets
                   strip.background = element_blank(),
                   strip.text = element_text(size = 6 ,face="bold"),
                   #background
                   panel.background = element_rect(fill = "lightblue")
)


# Build inset
g1_inset <-  ggplotGrob(
  ggplot() +
    # Plot land
    geom_sf(data=foreign, fill="grey80", color="white", lwd=0.3) +
    geom_sf(data=usa, fill="grey80", color="white", lwd=0.3) +
    # Plot box
    annotate("rect", xmin=-122.6, xmax=-121, ymin=36.2, ymax=37.1, color="black", fill=NA, lwd=0.8) +
    # Label regions
    #geom_text(data=region_labels, mapping=aes(y=lat_dd, label=region), x= -124.4, hjust=0, size=2) +
    # Labels
    labs(x="", y="") +
    # Crop
    coord_sf(xlim = c(-124.5, -117), ylim = c(32.5, 42)) +
    # Theme
    theme_bw() + my_theme +
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


###MUST RUN FigS1_kelp_urchin_timeseries
g2_inset <- ggplotGrob(g1)

# Create the "Monterey" text label
monterey_label <- data.frame(
  x = c(-121.9, -122.01), # x-coordinate for the upper right corner
  y = c(36.657, 36.54),  # y-coordinate for the upper right corner
  label = c("Monterey \nBay", "Carmel \nBay")
)

 # Plot the county boundaries and site locations with non-overlapping labels
map <- ggplot() +
  geom_sf(data = ca_counties_mpen) +
  geom_sf(data = site_locations_sf) +
  ggrepel::geom_label_repel(
    data = site_locations_labels,
    aes(x = longitude, y = latitude, label = site),
    box.padding = 0.3,
    point.padding = 0.5,
    force = 18,
    size = 2,
    min.segment.length = 0.1,
    segment.color = "grey50"
  ) +
  # Plot cities
  geom_text(data=landmarks, mapping=aes(x=long, y=lat, label=place),
            size=4, fontface= "bold") +
  coord_sf(xlim = c(-122.08, -121.85), ylim = c(36.5, 36.66)) +
  ggsn::north(x.min = -122.05, x.max = -121.86, 
              y.min = 36.5, y.max = 36.66,
              location = "topright", 
              scale = 0.05, 
              symbol = 10)+
  theme_bw() + my_theme +
  xlab("Longitude")+
  ylab("Latitude")+
  # Add CA inset
  annotation_custom(grob = g1_inset, 
                    xmin = -121.88, 
                    xmax = -121.827,
                    ymin = 36.63) +
  ggsn::scalebar(x.min = -122.088, x.max = -121.88, 
                 y.min = 36.5, y.max = 36.66,
                 #anchor=c(x=-124.7,y=41),
                 location="bottomleft",
                 dist = 2.5, dist_unit = "km",
                 transform=TRUE, 
                 model = "WGS84",
                 st.dist=0.02,
                 st.size=3,
                 border.size=.5,
                 height=.02
  )+
  ggsn::north(x.min = -122.05, x.max = -121.85, 
              y.min = 36.5, y.max = 36.665,
              location = "topright", 
              scale = 0.05, 
              symbol = 10)+
  # Add plot inset
  annotation_custom(grob = g2_inset, 
                    xmin = -122.093, 
                    xmax = -121.97,
                    ymin = 36.59) +
  #add labels
  geom_text(data=monterey_label, mapping=aes(x=x, y=y, label=label),
            size=4, fontface= "bold") +
  #add arrow
  annotate("segment", x = -122, y = 36.54, xend = -121.96, yend = 36.54,
           arrow = arrow(type = "closed", length = unit(0.02, "npc")))

map


ggsave(map, filename=file.path(figdir, "Fig1_site_map_newv2.png"), bg = "white",
       width=7.5, height=6, units="in", dpi=600) 



