

#Joshua G. Smith; jossmith@mbayaq.org

rm(list=ls())


######
#required packages
librarian::shelf(tidyverse, sf, raster, shiny, tmap, terra, tidyterra, ggnewscale)

#set directories 
basedir <- "/Volumes/seaotterdb$/kelp_recovery/data"
figdir <- here::here("analyses","2incipient_forests","figures")
output <- here::here("analyses","2incipient_forests","output")


#read state
ca_counties_orig <- st_read(file.path(basedir, "gis_data/raw/ca_county_boundaries/s7vc7n.shp")) 

#read landsat raw
landsat_orig <- st_read(file.path(basedir, "kelp_landsat/processed/monterey_peninsula/landsat_mpen_1984_2023_points_withNAs.shp"))

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
#Step 2 - define max kelp extent for each focal region

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
#Step 4 - plot

# Theme
base_theme <-  theme(axis.text=element_text(size=4, color = "black"),
                     axis.title=element_text(size=6,color = "black"),
                     plot.tag=element_text(size=7,color = "black"),
                     plot.title=element_text(size=6,color = "black", face = "bold"),
                     # Gridlines
                     panel.grid.major = element_blank(), 
                     panel.grid.minor = element_blank(),
                     panel.background = element_blank(), 
                     axis.line = element_line(colour = "black"),
                     # Legend
                     legend.key.size = unit(0.3, "cm"), 
                     legend.key = element_rect(fill=alpha('blue', 0)),
                     legend.spacing.y = unit(0.1, "cm"),  
                     legend.text=element_text(size=4,color = "black"),
                     legend.title=element_text(size=4,color = "black"),
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
  labs(fill = expression("Kelp biomass" ~(kg~"/"~"30m"^{-2})))+
  #add land
  geom_sf(data = ca_counties_orig, fill = "gray", color = "gray80") +
  coord_sf(xlim = c(-121.99, -121.88), ylim = c(36.519, 36.645), crs = 4326)+
  labs(title = "2022 Q3 kelp cover")+
  #reorder legend
  #guides(
  #  fill = guide_legend(order = 2),  
  #  fill_new_scale = guide_legend(order = 3),  
  #  fill_new_scale1 = guide_legend(order = 1)  
  #) +  
  #add landmarks
  #geom_text(data=monterey_label, mapping=aes(x=x, y=y, label=label),
  #         size=3, fontface= "bold") +
  theme_bw() + base_theme +theme(#axis.text.y = element_blank(),
    plot.tag.position = c(-0.03, 1),
    axis.title=element_blank(),
    panel.background = element_rect(fill = "navyblue"))
p1


# Save the final figure
#ggsave(p1, filename = file.path(figdir, "FigX_incipient_forests.png"), 
 #      width = 4, height = 4, units = "in", dpi = 600)