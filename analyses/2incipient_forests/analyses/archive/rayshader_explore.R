

rm(list=ls())

library(rayshader)
library(raster)

# Load the GeoTIFF file
elevation_raster <- raster("/Users/jossmith/Documents/Data/land_data/2018_ca_fema_region9_Job971776/2018_ca_fema_region9_Job971776.tif")
#elevation_raster <- raster("/Users/jossmith/Documents/Data/land_data/USGS_13_n37w122_20240207.tif")

# Reduce the resolution of the raster data (e.g., factor of 4 for demonstration)
elevation_raster_agg <- aggregate(elevation_raster, fact=12) 

# Convert the aggregated elevation data to a matrix
elevation_matrix <- raster_to_matrix(elevation_raster_agg)
#plot(elevation_raster_agg)


# Assuming elevation_matrix has been defined as shown above
elevation_matrix%>%
sphere_shade(texture = "imhof1") %>%
rayshader::plot_3d(heightmap = elevation_matrix, 
                   zscale = 5, 
                   fov = 1, 
                   zoom = 0.7,
                   theta = 210, 
                   phi = 25, 
                   windowsize = c(800, 800), 
                   water = FALSE,
                   waterdepth = 3,
                   soil_levels = 16,
                   background = "#130669")

















