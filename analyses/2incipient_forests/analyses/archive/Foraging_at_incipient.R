
#Joshua G. Smith; jossmith@mbayaq.org

rm(list=ls())


######
#required packages
librarian::shelf(tidyverse, sf, raster, shiny, tmap, gganimate)

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

#read foraging data
forage_orig <- read_csv(file.path(basedir,"/foraging_data/processed/foraging_data_2016_2023.csv"))


################################################################################
#set incipient forests
plot_data <- final_data %>%
  #set recovery category
  mutate(Incipient = ifelse(site_name == "Cannery row"|
                              site_name == "Pescadero upcoast"|
                              site_name == "Carmel pinnacles" |
                              site_name == "Pescadero downcoast" |
                              site_name == "Stillwater cove"|
                              site_name =="Carmel river" |
                              site_name == "Monastery" |
                              site_name == "Whaler's cove","Yes","No"),
         Incipient = factor(Incipient, levels = c("Yes","No"))) %>%
  #drop sandy sites
  mutate(site_name = factor(site_name))

################################################################################
#Step 3 -- prep foraging data

forage_build1 <- forage_orig %>%
  filter (year >= 2016) %>%
  mutate(quarter = ifelse(month == 1 | month == 2 | month == 3, 1,
                          ifelse(month == 4 | month == 5 | month == 6,2,
                                 ifelse(month == 7 | month == 8 | month == 9,3,4))))  %>%
#filter urchin prey only
 filter(prey == "pur" | prey == "urc") 

focal_patch <-  forage_orig %>%
  filter (year >=2016)%>%
  mutate(quarter = ifelse(month == 1:3,1,
                          ifelse(month == 4:6,2,
                                 ifelse(month==7:9,3,4)))) %>%
  #filter urchin prey only
  # filter(prey == "pur" | prey == "mus") %>%
  #define focal patch
  group_by(bout) %>%
  summarize(n_dives = n()) %>%
  mutate(focal_patch = ifelse(n_dives > 3, "yes","no"))


#join

forage_build2 <- left_join(forage_build1, focal_patch, by="bout")


#aggregate data by bout

forage_build3 <- forage_build2 %>%
  group_by(year, quarter, bout, focal_patch, prey) %>%
  summarize(n_prey = sum(preynum),
            lat = mean(lat),
            long = mean(long)
  ) %>%
  filter(!(is.na(lat)))%>%
  #make spatial
  st_as_sf(.,coords = c("long","lat"), crs=4326) 


forage_plot_dat <- forage_build3 

st_crs(forage_plot_dat)

################################################################################
#plot clusters

#create grid

landsat_hclust <-plot_data %>% filter(year == 2023 & quarter == 3) %>%
  #transform to teale
  st_transform(crs=3310)

r <- terra::rast(landsat_hclust, res=30)  # Build a raster
vr_clust <- terra::rasterize(landsat_hclust, r,"Incipient") 


# Create the "Monterey" text label
monterey_label <- data.frame(
  x = c(-121.9, -121.98), # x-coordinate for the upper right corner
  y = c(36.64, 36.54),  # y-coordinate for the upper right corner
  label = c("Monterey \nBay", "Carmel \nBay")
)



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
#g1_inset


p1 <- ggplot() +
  geom_sf(data = ca_counties, fill = "gray", color = "gray80") +
  # Add clusters
  tidyterra::geom_spatraster(data = as.factor(vr_clust), na.rm = TRUE) +
  scale_fill_manual(values = c("forestgreen","lightblue"), na.value = "transparent") +
  #add foraging point
  geom_sf(data = forage_plot_dat, #aes(color = factor(year)),
          color = "black", fill ="black",
          size = 2) +
  labs(title = "", tag = "") +
  # Add landmarks
  #geom_text(data = monterey_label, mapping = aes(x = x, y = y, label = label),
  #        size = 3, fontface = "bold") +
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
  guides(fill = "none") +base_theme +
  coord_sf(xlim = c(-121.99, -121.88), ylim = c(36.519, 36.645), crs = 4326, datum=NA) 
  

p1


################################################################################
#build animation

##Create animation with points showing up one by one
plot_anim <- p1 +
  transition_states(states = factor(year), state_length = 0, wrap = FALSE) +
  enter_recolor(fill = "#f0f5f9") +
  shadow_mark(past = TRUE, alpha = 1, fill = "#3a6589")
plot_anim


##Render animation
animate(plot_anim, end_pause = 60,
        height = 200, width = 400) # a higher res img would not upload here :(




