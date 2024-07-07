

#Joshua G. Smith; jossmith@mbayaq.org

rm(list=ls())

################################################################################

# This code creates distance-based hierarchical clusters using landsat data. 
# The landsat data are imported as a points file, where each point is the center 
# of a 30x30 meter pixel. 

#The steps are:
#1: Import the landsat data and find the maximum kelp area for each point. 
#2: Perform hierarchical clustering with a 500m threshold. 
#3: Map the cluster assignments back to each point in the original landsat data. 
#4: Export cluster assignments. 

################################################################################

#Notes:
#1. need to check site naming to make sure set.seed worked


######
#required packages
librarian::shelf(tidyverse, sf, raster, shiny, tmap, terra, tidyterra, RColorBrewer)

#set directories 
basedir <- "/Volumes/seaotterdb$/kelp_recovery/data"
localdir <- "/Users/jossmith/Documents/Data/landsat"

#read state
ca_counties_orig <- st_read(file.path(basedir, "gis_data/raw/ca_county_boundaries/s7vc7n.shp")) 

#read landsat raw
landsat_orig <- st_read(file.path(localdir, "/processed/monterey_peninsula/landsat_mpen_1984_2023_points_withNAs.shp"))

# Get land
usa <- rnaturalearth::ne_states(country="United States of America", returnclass = "sf")
foreign <- rnaturalearth::ne_countries(country=c("Canada", "Mexico"), returnclass = "sf")

#read bathy
#bathy_10m <- st_read(file.path(basedir, "gis_data/raw/bathymetry/contours_5m/contours_5m.shp")) %>% filter(CONTOUR == "-10")



################################################################################
################################################################################
#Step 1 - define max kelp extent

#define kelp extent for 10 years before the MHW
na_dat <- landsat_orig %>% 
          filter(year > 2003 & year < 2014)%>%
            filter (latitude >= 36.510140 &
                                     latitude <= 36.670574, 
                                   area > 0)

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
height_threshold <- 500  # Adjust this based on cluster size
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
  filter(year == 2023 & quarter == 3) %>%
  #transform to teale
  st_transform(plot_dat, crs=3310)


################################################################################
#save

landsat_cluster_ID <- landsat_clust_extent %>% st_drop_geometry()

saveRDS(landsat_cluster_ID, file = file.path(basedir,"/kelp_landsat/processed/landsat_cluster_ID_123.Rds")) #last write 16 Feb 2024

