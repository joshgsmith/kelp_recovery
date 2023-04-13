rm(list=ls())

#required packages
librarian::shelf(tidyverse, sf, ggplot2, leaflet, tmap)

#set directories 
basedir <- "/Volumes/seaotterdb$/RESEARCH FORMS and FILES/Site Selection"

#read polygon dat
poly_dat <- st_read(file.path(basedir, "ATOS_allCATEN.shp")) %>% st_transform(crs=3310)

#read point dat
buoy_dat <- st_read(file.path(basedir, "stationsshapefile.shp"))%>% st_transform(crs=3310)

#read point dat
buoy_tab <- read.csv(file.path(basedir, "summerbouydata.csv"))
buoy_tab[buoy_tab == 200.00] <- NA
buoy_tab[buoy_tab == -Inf] <- NA

#read state
ca_counties <- st_read("/Volumes/seaotterdb$/kelp_recovery/data/gis_data/raw/ca_county_boundaries/s7vc7n.shp")


#################################################################################
#check data structure
str(poly_dat)
str(buoy_dat)
str(buoy_tab)

#################################################################################
#process polygon data

#calculate centroid

poly_cen <- poly_dat %>% st_centroid()

#################################################################################
#process buoy tab data

#make spatial

buoy_tab_spat <- buoy_dat %>% st_as_sf(crs=3310)


#################################################################################
#calculate distance

#this returns row names
out <- st_nearest_feature(poly_cen, buoy_dat)


#this returns the element-wise distances 
dist <- st_distance(poly_cen, buoy_tab_spat[out,], by_element=TRUE) %>%
          data.frame() %>% rename("dist"=".")

#join
buoy_out_build1 <- cbind(poly_cen, dist) %>%
              mutate(dist_m = as.numeric(word(dist,1)))

buoy_out_build2 <- cbind(buoy_out_build1, buoy_tab_spat[out,])


#################################################################################
#join distance with envr vars

#make sure STATNUM is factor
buoy_out_build2 <- buoy_out_build2 %>% mutate(STATNUM = factor(STATNUM))
buoy_tab <- buoy_tab %>% mutate(STATNUM = factor(STATNUM)) %>%
                #drop duplicates with join
                dplyr::select(!(c("STATNAM","LAT","LONG")))


buoy_out_build3 <- left_join(buoy_out_build2, buoy_tab, by=c("STATNUM")) %>% 
                    dplyr::select(!(c(geometry,geometry.1)))



#rejoni polygon area
poly_area <- poly_dat %>% dplyr::select(c("OBJECTID","geometry"))
buoy_out_build4 <- st_join(poly_area, buoy_out_build3, by="OBJECTID" ) %>%
                      dplyr::select(!(c(OBJECTID.y, dist))) %>%
                      rename("OBJECTID" = "OBJECTID.x")


#st_write(".shp")


#################################################################################
#static plot

# Setup theme
my_theme <-  theme(axis.text=element_text(size=6),
                   axis.title=element_text(size=8),
                   legend.text=element_text(size=6),
                   legend.title=element_text(size=8),
                   strip.text=element_text(size=8),
                   plot.title=element_text(size=10),
                   # Gridlines
                   panel.border = element_blank(),
                   panel.grid.major = element_blank(), 
                   panel.grid.minor = element_blank(),
                   panel.background = element_blank(), 
                   axis.line = element_line(colour = "black"),
                   # Legend
                   legend.background = element_rect(fill=alpha('blue', 0)))

buoy_out_build4 %>% 
  #filter focal region so that legend scale is appropriate
  filter(LAT >= 36.50 & LAT <= 37.2) %>%
  ggplot()+
  #define envr var to fill polygons
  geom_sf(aes(fill=wvhtmax), col="black", show.legend=T)+
  scale_fill_gradient(low="blue", high="indianred",na.value = "gray")+
  #add land
  geom_sf(data=ca_counties, fill="tan",  lwd=0.3) +
  coord_sf(xlim = c(-122.3, -121.75), ylim = c(36.50, 37.2), crs=4326)+
  my_theme

#################################################################################
#testing dynamic plot

#select focal region for testing (can make statwide extent by commenting out the filter)
dyn_poly <- buoy_out_build4  %>% 
  #filter focal region so that legend scale is appropriate
  filter(LAT >= 36.50 & LAT <= 37.2)

tm_polygons <- tm_shape(dyn_poly) +  tm_fill("wvhtmax") + tm_polygons() 


##generate interactive plot. Click on cell to view value
tmap_leaflet(tm_polygons)





