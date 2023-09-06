

#Joshua G. Smith; jossmith@mbayaq.org

rm(list=ls())

####

#Notes:
#1. need to check site naming to make sure set.seed worked


######
#required packages
librarian::shelf(tidyverse, sf, raster, shiny, tmap, terra, tidyterra, RColorBrewer)

#set directories 
basedir <- "/Volumes/seaotterdb$/kelp_recovery/data"

#read state
ca_counties_orig <- st_read(file.path(basedir, "gis_data/raw/ca_county_boundaries/s7vc7n.shp")) 

#read landsat raw
landsat_orig <- st_read(file.path(basedir, "kelp_landsat/processed/monterey_peninsula/landsat_mpen_1984_2023_points_withNAs.shp"))


# Get land
usa <- rnaturalearth::ne_states(country="United States of America", returnclass = "sf")
foreign <- rnaturalearth::ne_countries(country=c("Canada", "Mexico"), returnclass = "sf")


################################################################################
################################################################################
#Step 1 - define max kelp extent

#define 0 cover as historical kelp footprint
na_dat <- landsat_orig %>% filter (latitude >= 36.510140 &
                                     latitude <= 36.670574, 
                                   biomass == 0)

################################################################################
#Step 2 - determine clustering

#find distinct landsat points across all years
clust_dat <- na_dat %>%
              distinct(latitude, longitude, geometry) 

# Create a spatial distance matrix
dist_matrix <- st_distance(clust_dat)

# Perform hierarchical clustering
set.seed(1985)
hierarchical_clusters <- hclust(as.dist(dist_matrix))

# Cut the tree to define clusters (adjust the height threshold as needed)
height_threshold <- 1500  # Adjust this based on cluster size
cluster_assignment <- cutree(hierarchical_clusters, h = height_threshold)

# Add cluster assignments back to the original 'sf' object
clust_dat$cluster <- cluster_assignment

# Visualize the clustered polygons
plot(clust_dat, col = clust_dat$cluster)

# Calculate the number of clusters assigned
length(unique(cluster_assignment))


################################################################################
#Step 3 - map clusters back to original landat data

landsat_cluster <- st_join(landsat_orig, clust_dat) 
                    
#clean up
landsat_cluster_clean <- landsat_cluster %>%
                          rename(latitude = latitude.x,
                                 longitude = longitude.x)%>%
                          dplyr::select(-latitude.y, -longitude.y)

#test if it worked for focal year and quarter

landsat_clust_extent <- landsat_cluster_clean %>% 
                      #select monterey area
                      filter(latitude >= 36.510140 &
                                latitude <= 36.670574) %>%
                      filter(year == 2022 & quarter == 3) %>%
                      #transform to teale
                      st_transform(plot_dat, crs=3310)


################################################################################
#rename clusters with something useful

#determine center coords (THIS WILL BE IMPORTANT TO SAVE LATER)
cluster_position <- landsat_clust_extent %>%
                      st_drop_geometry()%>%
                      group_by(cluster)%>%
                      dplyr::summarize(lat = mean(latitude),
                                       long = mean(longitude))%>%
  ungroup() %>%
  st_as_sf(coords = c("long", "lat"), crs = 4326)

ggplot(data = cluster_position) +
  geom_sf() +  # Plot the spatial data
  geom_sf(data = ca_counties_orig, fill = "gray", color = "gray80") +
  geom_sf_text(aes(label = cluster), size = 8, color = "black")+
  coord_sf(xlim = c(-121.99, -121.88), ylim = c(36.519, 36.645), crs = 4326) 


#rename sites and create table
site_cluster_table <- landsat_clust_extent %>%
              mutate(site_name = case_when(
                cluster == 26 ~ "Cannery row",
                cluster == 25 ~ "Hopkins upcoast",
                cluster == 24 ~ "Hopkins downcoast",
                cluster == 23 ~ "Siren",
                cluster == 21 ~"Piños upcoast",
                cluster == 16 ~ "Piños",
                cluster == 15 ~ "Piños downcoast",
                cluster == 11 ~"Asilomar",
                cluster == 8 ~ "Spanish bay",
                cluster == 7 ~ "Point Joe",
                cluster == 5 ~ "China rock",
                cluster == 3 ~ "Fanshell",
                cluster == 1 ~ "Cypress",
                cluster == 2 ~ "Sunset",
                cluster == 4 ~ "Pescadero upcoast",
                cluster == 6 ~ "Carmel pinnacles",
                cluster == 9 ~ "Pescadero downcoast",
                cluster == 14 ~ "Stillwater cove",
                cluster == 18 ~ "Arrowhead point",
                cluster == 19 ~ "Carmel main",
                cluster == 20 ~ "Carmel river",
                cluster == 22 ~ "Monastery",
                cluster == 17 ~ "Whaler's cove",
                cluster == 12 ~ "Sea lion point",
                cluster == 10 ~ "South lobos",
                TRUE ~ as.character(cluster)
              ),
              #set cluster order
              site_num = case_when(
              cluster == 26 ~ 1,
              cluster == 25 ~ 2,
              cluster == 24 ~ 3,
              cluster == 23 ~ 4,
              cluster == 21 ~ 5,
              cluster == 16 ~ 6,
              cluster == 15 ~ 7,
              cluster == 11 ~ 8, 
              cluster == 8 ~ 9,
              cluster == 7 ~ 10,
              cluster == 5 ~ 11,
              cluster == 3 ~ 12,
              cluster == 1 ~ 13,
              cluster == 2 ~ 14,
              cluster == 4 ~ 15,
              cluster == 6 ~ 16,
              cluster == 9 ~ 17,
              cluster == 14 ~ 18,
              cluster == 18 ~ 19,
              cluster == 19 ~ 20,
              cluster == 20 ~ 21,
              cluster == 22 ~ 22,
              cluster == 17 ~ 23,
              cluster == 12 ~ 24,
              cluster == 10 ~ 25,
              TRUE ~ NA)
              ) %>%
            #filter to data extent
            filter(!(is.na(site_num)))%>%
            #select vars
            dplyr::select(site_name, site_num, latitude, longitude) # each row is a distinct landsat point


################################################################################
#save

landsat_cluster_ID <- site_cluster_table %>% st_drop_geometry()

saveRDS(landsat_cluster_ID, file = file.path(basedir,"/kelp_landsat/processed/landsat_cluster_ID.Rds"))

################################################################################
#rasterize

#create grid
r <- rast(landsat_clust_extent, res=30)  # Builds a blank raster of given dimensions and resolution  
vr <- rasterize(landsat_clust_extent, r,"cluster", resolution = 30) 


################################################################################
#

# Theme
base_theme <-  theme(axis.text=element_text(size=6, color = "black"),
                     axis.title=element_text(size=7,color = "black"),
                     plot.tag=element_text(size=7,color = "black"),
                     plot.title=element_text(size=7,color = "black", face = "bold"),
                     # Gridlines
                     panel.grid.major = element_blank(), 
                     panel.grid.minor = element_blank(),
                     panel.background = element_blank(), 
                     axis.line = element_line(colour = "black"),
                     # Legend
                    # legend.key.size = unit(0.3, "cm"), 
                    # legend.key = element_rect(fill=alpha('blue', 0)),
                    # legend.spacing.y = unit(0.1, "cm"),  
                     legend.text=element_text(size=5,color = "black"),
                     legend.title=element_text(size=6,color = "black"),
                     #legend.key.height = unit(0.1, "cm"),
                     legend.background = element_rect(fill=alpha('blue', 0)),
                     #facets
                     strip.text = element_text(size=6, face = "bold",color = "black", hjust=0),
                     strip.background = element_blank())



# Create the "Monterey" text label
monterey_label <- data.frame(
  x = c(-121.9, -121.98), # x-coordinate for the upper right corner
  y = c(36.64, 36.54),  # y-coordinate for the upper right corner
  label = c("Monterey \nBay", "Carmel \nBay")
)


# Choose a color palette with distinct colors
custom_palette <- brewer.pal(n = 8, name = "Dark2")

# Expand the palette to 28 colors by repeating it
custom_palette <- rep(custom_palette, length.out = 28)

# Set "NA" values to be transparent
custom_palette[is.na(custom_palette)] <- "transparent"

p1 <- ggplot() +
  # Add clusters
  tidyterra::geom_spatraster(data = as.factor(vr), na.rm = TRUE) +
  scale_fill_manual(values = custom_palette, na.value = "transparent") +
  geom_sf(data = ca_counties_orig, fill = "gray", color = "gray80") +
  coord_sf(xlim = c(-121.99, -121.88), ylim = c(36.519, 36.645), crs = 4326) +
  labs(title = "2022", tag = "B") +
  # Add landmarks
  geom_text(data = monterey_label, mapping = aes(x = x, y = y, label = label),
            size = 3, fontface = "bold") +
  theme_bw() + base_theme + theme(axis.text.y = element_blank(),
                                  plot.tag.position = c(-0.03, 1),
                                  axis.title = element_blank())+
  guides(fill="none")
p1



