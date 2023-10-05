

#Joshua G. Smith; jossmith@mbayaq.org

rm(list=ls())


######
#required packages

librarian::shelf(tidyverse, sf, raster, terra)

basedir <- "/Volumes/seaotterdb$/kelp_recovery/data"
figdir <- here::here("analyses","5community_regulation","figures")

#read landsat dat
landsat_orig <- st_read(file.path(basedir,"kelp_landsat/processed/monterey_peninsula/landsat_mpen_1984_2023_points_withNAs.shp"))

#read foraging data
forage_orig <- read_csv(file.path(basedir,"/foraging_data/processed/foraging_data_2016_2023.csv"))

#read bathy
bathy_5m <- st_read(file.path(basedir, "gis_data/raw/bathymetry/contours_5m/contours_5m.shp")) %>% filter(CONTOUR == "-5" | CONTOUR == "-10")

#read state
ca_counties_orig <- st_read(file.path(basedir, "gis_data/raw/ca_county_boundaries/s7vc7n.shp")) 

# Get land
usa <- rnaturalearth::ne_states(country="United States of America", returnclass = "sf")
foreign <- rnaturalearth::ne_countries(country=c("Canada", "Mexico"), returnclass = "sf")


################################################################################
#determine max kelp extent as proxy for shallow subtidal


# transform landsat data to Teale Albers
rast_build1 <- st_transform(landsat_orig, crs = 3310) %>% 
  mutate(biomass = ifelse(biomass == 0, NA, biomass))  %>% filter (year == 2022) 


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
#Step 2 - extract mussel dives

forage_build1 <- forage_orig %>% 
                    #filter mussel dives
                    filter(prey == "mus") %>%
                    #select bout locations
                    dplyr::select(year, month, day, bout, lat, long) %>% distinct() %>%
                    filter(!is.na(lat))%>%
                    st_as_sf(coords = c("long","lat"), crs = 4326)


################################################################################
#Step 3 - plot


base_theme <-  theme(axis.text=element_text(size=7, color = "black"),
                     axis.text.x = element_blank(),
                     axis.text.y = element_blank(),
                     axis.title=element_text(size=8,color = "black"),
                     legend.text=element_text(size=7,color = "black"),
                     legend.title=element_text(size=8,color = "black"),
                     plot.tag=element_text(size=10,color = "black"),
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

# Create the "Monterey" text label
monterey_label <- data.frame(
  x = c(-121.9, -121.97), # x-coordinate for the upper right corner
  y = c(36.64, 36.54),  # y-coordinate for the upper right corner
  label = c("Monterey \nBay", "Carmel \nBay")
)




p1 <- ggplot() +
  # Add landmarks
  geom_text(data = monterey_label, mapping = aes(x = x, y = y, label = label),
            size = 3, fontface = "bold") +
  # Add CA inset
  annotation_custom(grob = g1_inset, 
                    xmin = -122.01, 
                    xmax = -121.96,
                    ymin = 36.625) +
  #add historic kelp extent
  tidyterra::geom_spatraster(data = kelp_na, na.rm = TRUE) +
  scale_fill_gradient(low = alpha("forestgreen", alpha=0.6),
                      high = alpha("forestgreen", alpha=0.6),
                      na.value = NA,
                      labels = c("Max kelp\nextent"))+
  guides(fill = guide_legend(override.aes = list(size = 1),
                             label.theme = element_text(color = "gray"))) +
  #labs(fill = "Max kelp \nextent")+
  guides(
    fill = guide_legend(
      override.aes = list(size = 1),  # Adjust the legend point size as needed
      title = NULL  # Remove the legend title
    )
  )+
  #add foraging bouts for mussels
  ggnewscale::new_scale_fill()+
  geom_sf(
    data = forage_build1 %>% filter(year == 2017),
    aes(fill = "Mussel forage \nbout"),
    size=0.1
  ) +
  labs(
    fill = NULL,  # Remove the title from the legend
    guide = guide_legend(title = NULL)  # Remove the legend title
  )+
  guides(
    fill = guide_legend(
      override.aes = list(size = 3),  # Adjust the legend point size as needed
      title = NULL  # Remove the legend title
    )
  )+
  #add land
  geom_sf(data = ca_counties_orig, fill = "gray", color = "gray80") +
  geom_sf(data = bathy_5m, fill = "black", color = "black") +
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
  labs(title = "", tag = "B")+
  theme_bw() + base_theme + theme(#axis.text.y = element_blank(),
    #legend.position = "none",
   # plot.tag.position = c(0.05, 1), 
    axis.title=element_blank()) +  
  #adjust legend position
  theme(
    legend.position = "bottom",         # Place legend at the bottom
    legend.justification = "right",     # Align the legend to the right
    legend.box = "vertical",         
    legend.margin = margin(t = -60),   # Adjust the top margin to move the legend up
    plot.margin = margin(b = 90),         # Add extra margin at the bottom of the plot
    legend.box.spacing = unit(0.001, "cm")
  )

p1


#save
ggsave(p1, filename = file.path(figdir, "FigX_bout_locations2.png"), 
       width = 5, height = 4, units = "in", dpi = 600)

                    
  







