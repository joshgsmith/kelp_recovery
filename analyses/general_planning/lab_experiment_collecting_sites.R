#Joshua G. Smith; jossmith@mbayaq.org

rm(list=ls())

######
######
#required packages
librarian::shelf(tidyverse, sf, raster, shiny, tmap, terra, tidyterra, RColorBrewer)

#set directories 
basedir <- "/Volumes/seaotterdb$/kelp_recovery/data"
figdir <- here::here("analyses","cdfw_scp","figures")

#read state
ca_counties_orig <- st_read(file.path(basedir, "gis_data/raw/ca_county_boundaries/s7vc7n.shp")) 

#read landsat raw
landsat_orig <- st_read(file.path(basedir, "kelp_landsat/processed/monterey_peninsula/landsat_mpen_1984_2023_points_withNAs.shp"))

#read cluster area
landsat_hclust <- readRDS(file.path(basedir,"/kelp_landsat/processed/landsat_cluster_ID.Rds"))

#load MPAs
mpas_orig <- readRDS(file.path(basedir, "gis_data/processed/CA_MPA_polygons.Rds"))

# Get land
usa <- rnaturalearth::ne_states(country="United States of America", returnclass = "sf")
foreign <- rnaturalearth::ne_countries(country=c("Canada", "Mexico"), returnclass = "sf")

################################################################################
# Join landsat data with cluster ID

landsat_build1 <- left_join(landsat_orig, landsat_hclust, by = c("latitude","longitude")) %>%
  #drop sites outside of cluster assigments
  filter(!(is.na(site_name)))


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



rast_build3 <- st_transform(landsat_orig, crs = 3310) %>% 
  mutate(biomass = ifelse(biomass == 0, NA, biomass))  %>% filter (year == 2008 & quarter ==3) 
landsat_rast_2008 <- rasterize(rast_build3, r, field = "biomass", fun = mean)

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


#----------------------step 1 plot clusters-------------------------------------

#create grid

landsat_hclust <- landsat_build1 %>% filter(year == 2022 & quarter == 3) %>%
  #transform to teale
  st_transform(crs=3310)

r <- rast(landsat_hclust, res=30)  # Builds a blank raster of given dimensions and resolution  
vr_clust <- terra::rasterize(landsat_hclust, r,"site_name") 


# Create the "Monterey" text label
monterey_label <- data.frame(
  x = c(-121.9, -121.98), # x-coordinate for the upper right corner
  y = c(36.64, 36.54),  # y-coordinate for the upper right corner
  label = c("Monterey \nBay", "Carmel \nBay")
)
 
collect_site <- data.frame(
  longitude = c(-121.967079, -121.881841), # x-coordinate for the upper right corner
  latitude = c(36.567670, 36.603888),  # y-coordinate for the upper right corner
  site_name = c("Primary \nLone Cypress", "Secondary \nShale Reef")
)


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



p3 <- ggplot() +
  # Add CA inset
  annotation_custom(grob = g1_inset, 
                    xmin = -122.01, 
                    xmax = -121.96,
                    ymin = 36.625) +
  geom_sf(data=mpas_orig, fill=ifelse(mpas_orig$type == "SMR","indianred1","lightblue"),alpha=0.6, color="black") +
  #add historic kelp extent
  tidyterra::geom_spatraster(data = kelp_na, na.rm = TRUE) +
  scale_fill_gradient(low = alpha("#7286AC",0.4),
                      high = alpha("#7286AC",0.4),
                      na.value = NA)+
  guides(fill = guide_legend(override.aes = list(size = 3),
                             label.theme = element_text(color = "white"))) +
  labs(fill = "Max kelp \nextent")+
  #add observed landsat for 2022 Q3
  ggnewscale::new_scale_fill()+
  tidyterra::geom_spatraster(data = landsat_rast_2022, na.rm = TRUE) +
  scale_fill_gradient2(low = "navyblue",mid="#1B9C00",high = "#FFFF79",
                       na.value = NA)+
  #scale_fill_viridis_c(na.value = NA)+
  labs(fill = expression("Kelp biomass" ~(Kg~"/"~"30m"^{-2})))+
  #add land
  geom_sf(data = ca_counties_orig, fill = "gray", color = "gray80") +
  # Add landmarks
  geom_text(data = monterey_label, mapping = aes(x = x, y = y, label = label),
            size = 3, fontface = "bold") +
  # Add CA inset
  annotation_custom(grob = g1_inset, 
                    xmin = -122.01, 
                    xmax = -121.96,
                    ymin = 36.625) +
  #add sites
  geom_point(
    data = collect_site,
    aes(x = longitude, y = latitude),
    shape = 19,            # Use shape 19 for round points
    size = 2,              # Adjust point size if needed
    color = "black",      # Point color
  ) +
  # Add jittered labels for site_name with boxes
  ggrepel::geom_label_repel(
    data = collect_site,
    aes(x = longitude, y = latitude, label = site_name),
    #box.padding = 1,
    #point.padding = 20,
    force = 5,
    size = 2,
    nudge_x = 0.02,
    nudge_y = 0.01,
    min.segment.length = unit(0,'lines'),
    segment.color = "black",
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
  coord_sf(xlim = c(-121.99, -121.88), ylim = c(36.519, 36.645), crs = 4326)+
  labs(title = "2022 Q3")+
  theme_bw() + base_theme +theme(#axis.text.y = element_blank(),
    plot.tag.position = c(-0.03, 1),
    axis.title=element_blank(),
    panel.background = element_rect(fill = "white"))
p3




#save
ggsave(p3, filename = file.path(figdir, "collecting_sites_for_lab_exp.png"), 
       width = 5, height = 6, units = "in", dpi = 600)



