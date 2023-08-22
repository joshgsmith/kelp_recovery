#Joshua G. Smith, jossmith@mbayaq.org

rm(list=ls())

#required packages
librarian::shelf(tidyverse, tidync, sf, rnaturalearth, usethis, raster, tidyterra)

#set directories 
basedir <- "/Volumes/seaotterdb$/kelp_recovery/data"
figdir <- here::here("analyses","2forage_behavior","figures")

#read landsat dat
landsat_dat <- st_read(file.path(basedir, "kelp_landsat/processed/monterey_peninsula/landsat_mpen_1984_2022_points.shp"))

#read state
ca_counties <- st_read(file.path(basedir, "gis_data/raw/ca_county_boundaries/s7vc7n.shp")) 

#read foraging data
forage_dat <- read.csv(file.path(basedir, "foraging_data/processed/foraging_data_2016_2023.csv")) 

#read forage metadata
forage_meta <- read.csv(file.path(basedir, "foraging_data/processed/forage_metadata.csv")) 


################################################################################
#Step 1 -- process landsat data for join

mpen_kelp <- landsat_dat %>% 
              #filter years that match otter data
              filter(year >= 2016 &
                       #filter data for central coast only
                       latitude >= 36.510140 &
                       latitude <= 36.670574) %>%
  #drop null dat
  filter(!(biomass == "0" | is.na(biomass)))


#transform landsat data to Teale Albers
t_dat <- st_transform(mpen_kelp, crs=3310)

#create grid for raster
r <- terra::rast(t_dat, res=30)  # Builds a blank raster of given dimensions and resolution  

#rasterize
vr <- rasterize(t_dat, r,"biomass", resolution = 30) 

################################################################################
#Step 2 -- process foraging dat for join

forage_build1 <- forage_dat %>%
                  mutate(quarter = ifelse(month == 1 | month == 2 | month == 3, 1,
                                          ifelse(month == 4 | month == 5 | month == 6,2,
                                                 ifelse(month == 7 | month == 8 | month == 9,3,4)))) %>%
                  #filter urchin prey only
                  filter(prey == "pur") 
                  
focal_patch <-  forage_dat %>%
                  mutate(quarter = ifelse(month == 1:3,1,
                          ifelse(month == 4:6,2,
                                 ifelse(month==7:9,3,4)))) %>%
                    #filter urchin prey only
                  filter(prey == "pur") %>%
                  #define focal patch
                  group_by(bout) %>%
                  dplyr::summarize(n_dives = n()) %>%
                  dplyr::mutate(focal_patch = ifelse(n_dives > 3, "yes","no"))


#join

forage_build2 <- left_join(forage_build1, focal_patch, by="bout")



#aggregate data by bout

forage_build3 <- forage_build2 %>%
                    group_by(year, quarter, bout, focal_patch) %>%
                    #summarize(n_prey = sum(preynum),
                     #         lat = mean(lat),
                      #        long = mean(long)
                       #       ) %>%
                  filter(!(is.na(lat)))%>%
                  #make spatial
                  st_as_sf(., coords = c("long","lat"), crs=4326)


forage_build4 <- st_transform(forage_build3, crs=3310)


################################################################################
#Step 3 -- calculate distance to nearest point


# Calculate the distance matrix between forage_build4 and t_dat
dist_matrix <- st_distance(forage_build4, t_dat)

# Initialize an empty data frame to store the results
result <- data.frame()

# Initialize a vector to store the distances
distances <- c()

# Loop through each unique combination of 'year' and 'quarter' in forage_build4
for (unique_year in unique(forage_build4$year)) {
  for (unique_quarter in unique(forage_build4$quarter)) {
    # Filter for the specific 'year' and 'quarter' combination
    forage_build4_subset <- forage_build4[forage_build4$year == unique_year & forage_build4$quarter == unique_quarter, ]
    
    # Find the corresponding rows in t_dat for the same 'year' and 'quarter' combination
    t_dat_subset <- t_dat[t_dat$year == unique_year & t_dat$quarter == unique_quarter, ]
    
    # Calculate the distance matrix between the filtered forage_build4_subset and t_dat_subset
    dist_matrix_subset <- st_distance(forage_build4_subset, t_dat_subset)
    
    # Find the indices of the closest points for this subset
    closest_indices <- apply(dist_matrix_subset, 1, function(row) {
      min_dist <- min(row, na.rm = TRUE)
      closest_idx <- which(row == min_dist)
      closest_idx[1]  # If multiple points have the same closest distance, take the first one
    })
    
    # Create a new data frame with the closest points from t_dat_subset
    closest_points <- t_dat_subset[closest_indices, ]
    
    # Calculate the distances between the points
    distances_subset <- apply(dist_matrix_subset, 1, function(row) {
      min_dist <- min(row, na.rm = TRUE)
      min_dist
    })
    
    # Append the distances to the overall distances vector
    distances <- c(distances, distances_subset)
    
    # Merge the closest_points data frame with forage_build4_subset based on 'year' and 'quarter'
    result_subset <- st_join(forage_build4_subset, closest_points, join = st_nearest_feature)
    
    # Append the subset result to the overall result
    result <- rbind(result, result_subset)
  }
}

# Add the distances column to the result data frame
result$distance_to_closest <- distances


################################################################################
#Step 4 -- plot


plot_dat <- result %>% mutate(
  observed_within_50 = ifelse(canopy == 1 | canopy == 0.5,"yes","no")
)


my_theme <-  theme(axis.text=element_text(size=6, color = "black"),
                   axis.text.y = element_text(angle = 90, hjust = 0.5, color = "black"),
                   axis.title=element_text(size=8, color = "black"),
                   plot.tag=element_text(size= 8, color = "black"), #element_text(size=8),
                   plot.title =element_text(size=7, face="bold", color = "black"),
                   # Gridlines 
                   panel.grid.major = element_blank(), 
                   panel.grid.minor = element_blank(),
                   panel.background = element_blank(), 
                   axis.line = element_line(colour = "black"),
                   # Legend
                   legend.key = element_blank(),
                   legend.background = element_rect(fill=alpha('blue', 0)),
                   legend.key.height = unit(1, "lines"), 
                   legend.text = element_text(size = 6, color = "black"),
                   legend.title = element_text(size = 7, color = "black"),
                   #legend.spacing.y = unit(0.75, "cm"),
                   #facets
                   strip.background = element_blank(),
                   strip.text = element_text(size = 6 ,face="plain", hjust=0, color = "black"),
)

# Create a subset of the data with the "yes" and "no" points
yes_points <- subset(plot_dat, observed_within_50 == "yes")
no_points <- subset(plot_dat, observed_within_50 == "no")


# Calculate the counts and proportions for "yes" points
yes_counts <- table(yes_points$distance_to_closest <= 50)
yes_proportions <- yes_counts / sum(yes_counts)

# Calculate the counts and proportions for "no" points
no_counts <- table(no_points$distance_to_closest > 50)
no_proportions <- no_counts / sum(no_counts)

# Create a data frame to use for plotting
plot_data <- data.frame(
  canopy = c("Inside (aligned)", "Inside (out of sample)", "Outside (out of sample)", "Outside (aligned)"),
  count = c(yes_counts[TRUE], yes_counts[FALSE], no_counts[TRUE], no_counts[FALSE]),
  proportion = c(yes_proportions[TRUE], yes_proportions[FALSE], no_proportions[TRUE], no_proportions[FALSE])
)

g <- ggplot(plot_data, aes(x = canopy, y = count, fill = factor(canopy))) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = sprintf("%.2f", proportion)), position = position_stack(vjust = 0.5)) +
  labs(x = "Location", y = "No. observations") +
  scale_fill_manual(values = c("Inside (aligned)" = "#91bfdb", "Inside (out of sample)" = "#3182bd",
                               "Outside (out of sample)" = "#c51b7d", "Outside (aligned)" = "#f768a1")) +
  guides(fill = FALSE) +  # Remove the legend for the "fill" aesthetic
  theme_bw() + my_theme

g

ggsave(filename = file.path(figdir, "FigX_distance_to_kelp.png"), plot = g, 
       width = 7, height = 6, bg = "white", units = "in", dpi = 600)







