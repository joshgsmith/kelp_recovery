

#Joshua G. Smith; jossmith@mbayaq.org

rm(list=ls())


######
#required packages
librarian::shelf(tidyverse, sf,rgeos,sp)

#set directories 
basedir <- "/Volumes/seaotterdb$/kelp_recovery/data/census_data"
censusdir <- file.path(basedir,"annual_surveys/raw")
output <-file.path(basedir,"annual_surveys/processed")


################################################################################
#Read USGS 1985-2014 shapefiles

#------------------------step 1, uncompress shapefiles--------------------------
zip_dir <- file.path(basedir,"annual_surveys/raw/zip_files")
uncomp_dir <- file.path(basedir,"annual_surveys/raw/uncompressed_files")


#usgs_1985_2014
#zip_files <- list.files(file.path(zip_dir,"usgs_1985_2014"), pattern = "\\.zip$", full.names = TRUE)

#for (zip_file in zip_files) {
#  unzip(zip_file, exdir = uncomp_dir)
#}


#------------------------step 2, read shapefiles--------------------------------


# Function to identify unique columns in a list of data frames
getUniqueColumns <- function(data_list) {
  unique_cols <- unique(unlist(lapply(data_list, names)))
  return(unique_cols)
}

# Initialize a list to store unique columns
unique_columns <- NULL

# Initialize a list to store processed shapefiles
shapefile_list <- list()

# Iterate over subdirectories and process shapefiles
subdirectories <- list.dirs(uncomp_dir, recursive = FALSE)
for (subdir in subdirectories) {
  shapefiles <- list.files(subdir, pattern = "^Census_sum.*\\.shp$", full.names = TRUE)
  
  if (length(shapefiles) > 0) {
    subdir_year <- as.integer(substring(subdir, regexpr("\\d{4}", subdir), regexpr("\\d{4}", subdir) + 3))
    
    # Read each shapefile into a list with the specified CRS
    shapefile_data <- lapply(shapefiles, function(file) {
      sf_data <- sf::st_read(file)
      
      # Add a "Year_file" column with the extracted year
      sf_data$Year_file <- subdir_year
      
      return(sf_data)
    })
    
    # Store the list of shapefiles in shapefile_list
    shapefile_list[[subdir]] <- shapefile_data
    
    # Update the unique columns
    unique_columns <- union(unique_columns, getUniqueColumns(shapefile_data))
  }
}

# Iterate over the processed shapefiles and ensure they have all unique columns
for (subdir in subdirectories) {
  if (length(shapefile_list[[subdir]]) > 0) {
    for (i in seq_along(shapefile_list[[subdir]])) {
      missing_cols <- setdiff(unique_columns, names(shapefile_list[[subdir]][[i]]))
      if (length(missing_cols) > 0) {
        shapefile_list[[subdir]][[i]][missing_cols] <- NA
      }
    }
  }
}


# Bind rows of all shapefiles into a single shapefile
final_shapefile <- do.call(rbind, unlist(shapefile_list, recursive = FALSE))

#------------------------step 3, inspect and clean------------------------------

#inspect

na.year <- final_shapefile %>% filter(is.na(Year_file)) #should be 0
unmatched_year <- final_shapefile %>% filter(Year_file != Year) #should be 0
st_crs(final_shapefile)

#clean up
census_dat <- final_shapefile %>%
  rownames_to_column(var = "RowName") %>%
  dplyr::select(-RowName, -Year)%>%
  rename(Year = Year_file)%>%
  janitor::clean_names()

#export

st_write(census_dat, file.path(output, "usgs_1985_2019_annual_census.shp")) 
s



