#Joshua G. Smith; jossmith@mbayaq.org

rm(list=ls())

######
######
#required packages
librarian::shelf(tidyverse, sf, raster, shiny, tmap, terra, tidyterra, RColorBrewer)

#set directories 
basedir <- "/Volumes/seaotterdb$/kelp_recovery/data"
localdir <- "/Users/jossmith/Documents/Data/landsat/processed" #local directory
figdir <- here::here("analyses","2incipient_forests","figures")
output <- here::here("analyses","2incipient_forests","output")

#read state
ca_counties_orig <- st_read(file.path(basedir, "gis_data/raw/ca_county_boundaries/s7vc7n.shp")) 

#read landsat raw
landsat_orig <- st_read(file.path(localdir, "/monterey_peninsula/landsat_mpen_1984_2023_points_withNAs.shp"))


# Get land
usa <- rnaturalearth::ne_states(country="United States of America", returnclass = "sf")
foreign <- rnaturalearth::ne_countries(country=c("Canada", "Mexico"), returnclass = "sf")


################################################################################
#Step 1 - raster for focal year and quarter

# transform landsat data to Teale Albers
rast_build1 <- st_transform(landsat_orig, crs = 3310) %>% 
  mutate(area = ifelse(area == 0, NA, area))  %>% filter (year == 2023 & quarter ==3) 

#define blank raster
r <- rast(rast_build1, res=30)

landsat_rast_2023 <- rasterize(rast_build1, r, field = "area", fun = mean)


################################################################################
#Step 2 - define max kelp extent

#select MPEN as focal region 
plot_dat <- rast_build1 %>% filter (latitude >= 36.510140 &
                                      latitude <= 36.670574) %>% st_transform(crs=3310)

#define 0 cover as historical kelp footprint
na_dat <- landsat_orig %>% filter (latitude >= 36.510140 &
                                     latitude <= 36.670574, 
                                   area == 0)

#transform landsat data to Teale Albers
t_dat <- st_transform(plot_dat, crs=3310) %>% 
  mutate(area = ifelse(area==0,NA,area))

kelp_historic <- st_transform(na_dat, crs=3310) # %>% mutate(area = ifelse(area==0,1,NA))

#create grid
r <- rast(t_dat, res=30)  # Builds a blank raster of given dimensions and resolution  
vr <- rasterize(t_dat, r,"area", resolution = 30) 

kelp_na <- rasterize(kelp_historic, r,"area", resolution = 30) 


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
  tidyterra::geom_spatraster(data = landsat_rast_2023, na.rm = TRUE) +
  scale_fill_gradient2(low = "navyblue",mid="#1B9C00",high = "#FFFF79",
                       na.value = NA)+
  #add land
  geom_sf(data = ca_counties_orig, fill = "gray", color = "gray80") +
  #add landmarks
 # geom_text(data=monterey_label, mapping=aes(x=x, y=y, label=label),
  #          size=3, fontface= "bold") +
  coord_sf(xlim = c(-121.99, -121.88), ylim = c(36.519, 36.645), crs = 4326)+
  labs(title = "2023 Q3", tag = "")+
  theme_bw() + base_theme + theme(#axis.text.y = element_blank(),
    #legend.position = "none",
    plot.tag.position = c(0.05, 1), axis.title=element_blank())
p1

#write shapefiles for ArcGIS
library(raster)

# Write landsat_rast_2023 raster to a shapefile
#convert to vector format 
rast_my_spat <- raster(landsat_rast_2023)
poly_rast <- rasterToPolygons(rast_my_spat)
landsat_shape <- shapefile(poly_rast, file.path(basedir,"/kelp_landsat/processed/monterey_peninsula/ArcGIS_files/landsat_rast_2023.shp"))

# write as KML
landsat_shape <- st_read(file.path(basedir, "/kelp_landsat/processed/monterey_peninsula/ArcGIS_files/landsat_rast_2023.shp"))
st_write(landsat_shape, file.path(basedir, "/kelp_landsat/processed/monterey_peninsula/ArcGIS_files/landsat_rast_2023.kml"), driver = "KML")

dissolved_rast<- st_union(landsat_shape)
st_write(dissolved_rast, file.path(basedir, "/kelp_landsat/processed/monterey_peninsula/ArcGIS_files/kelp_rast_2023_polys.kml"), driver = "KML")




rast_my_spat <- raster(kelp_na)
poly_rast <- rasterToPolygons(rast_my_spat)
landsat_shape <- shapefile(poly_rast, file.path(basedir,"/kelp_landsat/processed/monterey_peninsula/ArcGIS_files/kelp_max_extent.shp"))

#write as KML
kelp_na <- st_read(file.path(basedir, "/kelp_landsat/processed/monterey_peninsula/ArcGIS_files/kelp_max_extent.shp"))
# Dissolve based on shared borders
dissolved_kelp <- st_union(kelp_na)

st_write(kelp_na, file.path(basedir, "/kelp_landsat/processed/monterey_peninsula/ArcGIS_files/kelp_max_extent.kml"), driver = "KML")
st_write(dissolved_kelp, file.path(basedir, "/kelp_landsat/processed/monterey_peninsula/ArcGIS_files/kelp_max_polys.kml"), driver = "KML")


#save
#ggsave(p1, filename = file.path(figdir, "area_2023_Q3.png"), 
 #      width = 5, height = 6, units = "in", dpi = 600)
