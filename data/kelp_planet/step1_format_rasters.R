




################################################################################
# About
# data processing script written by JG.Smith jossmith@mbayaq.org



################################################################################

rm(list=ls())

librarian::shelf(googledrive, tidyverse, raster, sf)
gs4_auth()

#set dir
datdir <- "/Volumes/seaotterdb$/kelp_recovery/data/MBA_kelp_forest_database/"

#set drive location
folder_url <- "https://drive.google.com/drive/folders/1MsMxxsPpCqexIGAemkDZEAIDwsHa7Z_e"

# List the files in folder
files <- drive_ls(as_id(folder_url))

# Filter for .tif files
tif_files <- files[grep("northCoast_refugia_\\d{4}\\.tif", files$name), ]

# Initialize an empty list to store rasters
raster_list <- list()

# Loop to read each .tif file into R's memory
for (i in 1:nrow(tif_files)) {
  # Extract year from the filename
  year <- gsub("northCoast_refugia_|\\.tif", "", tif_files$name[i])
  
  # Download to a temporary file and load directly as a raster
  temp_file <- tempfile(fileext = ".tif")
  drive_download(as_id(tif_files$id[i]), path = temp_file, overwrite = TRUE)
  raster_list[[year]] <- raster(temp_file)
}

# Example: Access and plot the raster for 2017
plot(raster_list[["2022"]], main = "North Coast Refugia 2022")


################################################################################
#Explore coordinate system projection before scaling

crs(raster_2022)
extent(raster_2022)
extent_to_test <- extent(425000, 450000, 4280000, 4320000) # Adjust to a smaller area within your range
raster_subset <- crop(raster_2022, extent_to_test)
raster_2022_wgs84_subset <- projectRaster(raster_subset, crs = sp::CRS(SRS_string = "EPSG:4326"))
plot(raster_2022_wgs84_subset)

#check res
res(raster_subset) #3 x 3 meters -- incredible!!

# Define the extent in UTM Zone 10N
xmin <- 404997
xmax <- 489285
ymin <- 4254090
ymax <- 4431240

# Create spatial points for the corners of the extent
utm_coords <- SpatialPoints(
  matrix(c(xmin, ymin, xmax, ymax), ncol = 2, byrow = TRUE),
  proj4string = CRS("+proj=utm +zone=10 +datum=WGS84 +units=m +no_defs")
)

# Transform to WGS84 (decimal degrees)
wgs84_coords <- spTransform(utm_coords, CRS("+proj=longlat +datum=WGS84"))

# View the transformed coordinates
coordinates(wgs84_coords)

#
plot(raster_2022_wgs84_subset, col = c("forestgreen", "white"), main = "")








