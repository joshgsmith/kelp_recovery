

#Joshua G. Smith; jossmith@mbayaq.org

rm(list=ls())


######
#required packages
librarian::shelf(tidyverse, sf, raster)

#set directories 
basedir <- "/Volumes/seaotterdb$/kelp_recovery/data"
figdir <- here::here("analyses","2incipient_forests","figures")
output <- here::here("analyses","2incipient_forests","output")

#read USGS Atos
atos_orig <- st_read(file.path(basedir, "gis_data/raw/ATOS/ATOS_polygon_teale83/ATOS_polygon_teale83.shp"))

#read state
ca_counties_orig <- st_read(file.path(basedir, "gis_data/raw/ca_county_boundaries/s7vc7n.shp")) 

# Get land
usa <- rnaturalearth::ne_states(country="United States of America", returnclass = "sf")
foreign <- rnaturalearth::ne_countries(country=c("Canada", "Mexico"), returnclass = "sf")

################################################################################
#inspect

st_crs(atos_orig)

################################################################################
#select Monterey Bay as focal region

atos_mon <- atos_orig %>%
  st_transform(st_crs(ca_counties_orig)) %>%
  # Crop to Monterey region
  st_crop(st_bbox(c(ymin = 36.519, xmin = -121.99, ymax = 36.7, xmax = -121.88)))

#inspect
plot(atos_mon)

ggplot() +
  geom_sf(data = atos_mon) +
  #add land
  geom_sf(data = ca_counties_orig, fill = "gray", color = "gray80") +
  coord_sf(xlim = c(-121.99, -121.88), ylim = c(36.519, 36.645), crs = 4326) +
  theme_bw()

################################################################################
#join 0-30 and 30-60 polygons



library(sf)


# Identify touching polygons across both DEPTH levels
touching_pairs <- st_touches(atos_mon$geometry)

# Extract unique touching pairs
unique_touching_pairs <- lapply(touching_pairs, function(ids) unique(atos_mon$ID[ids]))

# Combine touching pairs into lists of unique IDs
joined_list <- list()
for (i in seq_along(unique_touching_pairs)) {
  if (length(unique_touching_pairs[[i]]) > 1) {
    joined_geoms <- atos_mon %>%
      filter(ID %in% unique_touching_pairs[[i]]) %>%
      st_union()
    
    joined_geom_with_depth <- st_sf(
      geometry = joined_geoms,
      DEPTH = "Joined"
    )
    
    joined_list[[length(joined_list) + 1]] <- joined_geom_with_depth
  }
}

# Combine the joined polygons
joined_geoms <- bind_rows(joined_list, .id = "ID") %>%
  st_as_sf()

# Print the resulting joined geometries
print(joined_geoms)

plot(joined_geoms)












library(sf)

# Assuming atos_mon contains your cropped and transformed geometries with 'DEPTH' column

# Identify touching polygons across both DEPTH levels
touching_pairs <- st_touches(atos_mon$geometry)

# Combine touching pairs into lists of unique IDs
joined_list <- list()
processed_ids <- integer(0)  # To keep track of processed IDs

for (i in seq_along(touching_pairs)) {
  touching_ids <- unique(atos_mon$ID[touching_pairs[[i]]])
  
  # Check if any touching polygon is already processed
  if (any(touching_ids %in% processed_ids)) {
    next
  }
  
  # Extract DEPTH values of touching polygons
  touching_depths <- unique(atos_mon$DEPTH[touching_ids])
  
  # Check if DEPTH values are different
  if (length(touching_depths) > 1) {
    touching_geoms <- atos_mon %>%
      filter(ID %in% touching_ids) %>%
      st_make_valid()  # Ensure geometry validity
    
    joined_geoms <- st_union(touching_geoms$geometry)
    
    # Convert multipolygons to polygons
    if (st_geometry_type(joined_geoms) == "MULTIPOLYGON") {
      joined_geoms <- st_cast(joined_geoms, "POLYGON")
    }
    
    joined_geom_with_depth <- st_sf(
      geometry = joined_geoms,
      DEPTH = "Joined"
    )
    
    joined_list[[length(joined_list) + 1]] <- joined_geom_with_depth
    processed_ids <- c(processed_ids, touching_ids)
  }
}

# Combine the joined polygons
joined_geoms <- bind_rows(joined_list, .id = "ID") %>%
  st_as_sf()

# Print the resulting joined geometries
print(joined_geoms)


plot(joined_geoms)
