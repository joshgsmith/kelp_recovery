

#Joshua G. Smith; jossmith@mbayaq.org

rm(list=ls())


######
#required packages
librarian::shelf(tidyverse, sf, raster, shiny, tmap, terra, tidyterra, RColorBrewer)

#set directories 
basedir <- "/Volumes/seaotterdb$/kelp_recovery/data"
figdir <- here::here("analyses","2incipient_forests","figures")

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


#check number of clusters using elbow plot
threshold_values <- seq(500, 3500, by = 100)  
num_clusters <- numeric(length(threshold_values))

for (i in seq_along(threshold_values)) {
  cluster_assignment <- cutree(hierarchical_clusters, h = threshold_values[i])
  num_clusters[i] <- length(unique(cluster_assignment))
}

# Create a data frame
elbow_data <- data.frame(Threshold = threshold_values, NumClusters = num_clusters)

################################################################################
#plot

# Theme
base_theme <-  theme(axis.text=element_text(size=8, color = "black"),
                     axis.title=element_text(size=8,color = "black"),
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

# 
g <- ggplot(elbow_data, aes(x = Threshold, y = NumClusters)) +
  geom_line() +
  geom_point() +
  labs(x = "Height threshold", y = "Number of clusters") +
  geom_vline(xintercept = 1500, linetype = "dashed", color = "red") +
  annotate("text", x = 1550, y = max(num_clusters) - 1, 
           label = "Elbow point\nthreshold = 1800", 
           color = "red", hjust = 0, size=2) +
  geom_hline(yintercept = 28, linetype = "dashed", color = "black") +
  annotate("text", x = 200, y = 28 + 3, 
           label = "Cluster assignment = 26", 
           color = "black", hjust = 0, size=2) + theme_bw() + base_theme

g



ggsave(g, filename = file.path(figdir, "FigS1_elbow_plot.png"), 
       width = 4, height = 4, units = "in", dpi = 600)

