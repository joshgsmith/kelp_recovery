

#Joshua G. Smith; jossmith@mbayaq.org

rm(list=ls())


######
#required packages
librarian::shelf(tidyverse, sf, raster, shiny, tmap, terra, tidyterra, ggnewscale)

#set directories 
localdir <- "/Users/jossmith/Documents/Data/landsat/processed" #local directory
basedir <- "/Volumes/seaotterdb$/kelp_recovery/data" #remote directory
gisdir <- file.path(basedir, "gis_data/processed")
figdir <- here::here("analyses","2incipient_forests","figures")
output <- here::here("analyses","2incipient_forests","output")

#read MPAs
mpas <- readRDS(file.path(gisdir, "CA_MPA_polygons.Rds")) %>%
  dplyr::select(name, lat_dd, long_dd, type)

#read state
ca_counties_orig <- st_read(file.path(basedir, "gis_data/raw/ca_county_boundaries/s7vc7n.shp")) 

#read landsat raw
landsat_mpen <- st_read(file.path(basedir, "kelp_landsat/processed/monterey_peninsula/landsat_mpen_1984_2023_points_withNAs.shp"))
landsat_GFNMS <- st_read(file.path(localdir, "gfnms_sonoma/landsat_mpen_1984_2023_points_withNAs.shp"))
landsat_trin <- st_read(file.path(localdir, "trinidad/landsat_mpen_1984_2023_points_withNAs.shp"))

# Get land
usa <- rnaturalearth::ne_states(country="United States of America", returnclass = "sf")
foreign <- rnaturalearth::ne_countries(country=c("Canada", "Mexico"), returnclass = "sf")


################################################################################
#Step 1 - merge datasets

landsat_orig <- rbind(landsat_mpen, landsat_GFNMS, landsat_trin)


################################################################################
#Step 1 - raster for focal year and quarter

#2008
rast_build1 <- st_transform(landsat_orig, crs = 3310) %>% 
  mutate(biomass = ifelse(biomass == 0, NA, biomass))  %>% filter (year == 2008 & quarter ==3) 
#define blank raster
r <- rast(rast_build1, res=30)
landsat_rast_2008 <- rasterize(rast_build1, r, field = "area", fun = mean)

#2017
rast_build2 <- st_transform(landsat_orig, crs = 3310) %>% 
  mutate(biomass = ifelse(biomass == 0, NA, biomass))  %>% filter (year == 2017 & quarter ==3) 

landsat_rast_2017 <- rasterize(rast_build2, r, field = "area", fun = mean)


#2022
rast_build3 <- st_transform(landsat_orig, crs = 3310) %>% 
  mutate(biomass = ifelse(biomass == 0, NA, biomass))  %>% filter (year == 2022 & quarter ==3) 

landsat_rast_2022 <- rasterize(rast_build3, r, field = "area", fun = mean)



################################################################################
#Step 2 - define max kelp extent for each focal region


#2008
#select GF as focal region 
plot_dat_2008 <- rast_build1 %>% st_transform(crs=3310)

#define 0 cover as historical kelp footprint
na_dat <- landsat_orig %>% filter(area == 0)

#transform landsat data to Teale Albers
t_dat_2008 <- st_transform(plot_dat_2008, crs=3310) %>% 
  mutate(area = ifelse(area==0,NA,area))

kelp_historic <- st_transform(na_dat, crs=3310) # %>% mutate(biomass = ifelse(biomass==0,1,NA))

#create grid
r_2008 <- rast(t_dat_2008, res=30)  # Builds a blank raster of given dimensions and resolution  
vr_2008 <- rasterize(t_dat_2008, r_2008,"area", resolution = 30) 

kelp_na_2008 <- rasterize(kelp_historic, r_2008,"area", resolution = 30) 


#2017
#select GF as focal region 
plot_dat_2017 <- rast_build2 %>% st_transform(crs=3310)

#define 0 cover as historical kelp footprint
na_dat <- landsat_orig %>% filter(area == 0)

#transform landsat data to Teale Albers
t_dat_2008 <- st_transform(plot_dat_2008, crs=3310) %>% 
  mutate(area = ifelse(area==0,NA,area))

kelp_historic <- st_transform(na_dat, crs=3310) # %>% mutate(biomass = ifelse(biomass==0,1,NA))

#create grid
r_2008 <- rast(t_dat_2008, res=30)  # Builds a blank raster of given dimensions and resolution  
vr_2008 <- rasterize(t_dat_2008, r_2008,"area", resolution = 30) 

kelp_na_2008 <- rasterize(kelp_historic, r_2008,"area", resolution = 30) 




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


MPEN <- ggplot() +
  # Plot MPAs
  geom_sf(data=mpas, fill=ifelse(mpas$type == "SMR","indianred1","lightblue"), color="black", alpha=0.3) +
  #add historic kelp extent
  tidyterra::geom_spatraster(data = kelp_na_MPEN, na.rm = TRUE) +
  scale_fill_gradient(low = alpha("#8499C2",0.4),
                      high = alpha("#8499C2",0.4),
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
  labs(title = "2022")+
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
    panel.background = element_rect(fill = "white"))
MPEN



GFNMS <- ggplot() +
  # Plot MPAs
  geom_sf(data=mpas, fill=ifelse(mpas$type == "SMR","indianred1","lightblue"), color="black", alpha=0.3) +
  #add historic kelp extent
  tidyterra::geom_spatraster(data = kelp_na_GFNMS, na.rm = TRUE) +
  scale_fill_gradient(low = alpha("#8499C2",0.4),
                      high = alpha("#8499C2",0.4),
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
  coord_sf(xlim = c(-123.753617, -123.238646), ylim = c(38.508139, 38.953937), crs = 4326)+
  labs(title = "2022")+
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
    panel.background = element_rect(fill = "white"))
GFNMS








# Save the final figure
#ggsave(p1, filename = file.path(figdir, "FigX_incipient_forests.png"), 
 #      width = 4, height = 4, units = "in", dpi = 600)