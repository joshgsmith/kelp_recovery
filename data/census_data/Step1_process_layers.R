

#Joshua G. Smith; jossmith@mbayaq.org

rm(list=ls())


######
#required packages
librarian::shelf(tidyverse, sf, fs, janitor)

#set directories 
basedir <- "/Volumes/seaotterdb$/"
census_dir <- "RESEARCH FORMS AND FILES/Monterey monthly census/QuarterlyCensus/qc_raw"

################################################################################
#read data

# Get a list of all shapefile paths
shapefile_paths <- dir_ls(file.path(basedir, census_dir), recurse = TRUE, glob = "*.shp")

# Function to extract date from file name 
# some shapefiles do not have a date as an attribute, so we'll use the 
# date listed in the file name as add it as 'file_dt'
# Function to extract date from file path
extract_date <- function(path) {
  date_str <- sub('.*(\\d{8}).*', '\\1', path)
  as.Date(date_str, format = "%m%d%Y")
}

#inspect function to make sure it worked
extracted_dates <- lapply(shapefile_paths, extract_date)
print(extracted_dates)

# Read shapefiles into a list of sf objects and add 'file_dt' column
shapefile_list <- lapply(shapefile_paths, function(path) {
  sf_object <- st_read(path) %>% clean_names()
  sf_object$file_dt <- extract_date(path)
  sf_object
})

#inspect
View(shapefile_list[[14]])

#save geometry for later
save_crs <- st_crs(shapefile_list[[1]])

################################################################################
#standardize the names

# Extract column names from each sf object
list_of_column_names <- map(shapefile_list, colnames)

# Find the maximum number of columns among the lists
max_cols <- max(lengths(list_of_column_names))

# Pad the lists with missing column names using data.frame()
list_of_padded_column_names <- map(list_of_column_names, ~ data.frame(t(.x), stringsAsFactors = FALSE))

# Bind the rows to create the dataframe
column_names_df <- bind_rows(list_of_padded_column_names, .id = "ListName")

View(column_names_df)

# Drop geometry column from each sf object
shapefile_list <- map(shapefile_list, ~ .x %>% st_drop_geometry())

# Mapping of old column names to new standardized names
column_name_mapping <- c(
  lon = "long",
  x = "long",
  lat = "lat",
  y = "lat",
  numotters = "num_indepen",
  numindepend = "num_indepen",
  f_indepen = "num_indepen",
  numlgpups = "num_pups_lar",
  f_large_p = "num_pups_lar",
  numsmpups = "num_pups_sma",
  f_small_p = "num_pups_sma"
)

# Function to standardize column names
standardize_names <- function(sf_object) {
  renamed_sf <- sf_object %>%
    rename_with(~ ifelse(.x %in% names(column_name_mapping), column_name_mapping[.x], .x),
                cols = all_of(names(column_name_mapping)))
  return(renamed_sf)
}


# Apply the function to each sf object in the list
shapefile_list2 <- map(shapefile_list, ~ standardize_names(.x))

#recheck names
list_of_column_names2 <- map(shapefile_list2, colnames)
max_cols2 <- max(lengths(list_of_column_names2))
list_of_padded_column_names2 <- map(list_of_column_names2, ~ data.frame(t(.x), stringsAsFactors = FALSE))
column_names_df2 <- bind_rows(list_of_padded_column_names2, .id = "ListName")

# Combine all the unique column names from shapefile_list2
all_column_names <- unique(unlist(map(shapefile_list2, colnames)))




################################################################################
#merge into single dataframe

# Create a function to ensure consistent data types
convert_to_character <- function(df) {
  df <- df %>%
    mutate(across(everything(), as.character))
  return(df)
}

# Convert all dataframes to character data type
shapefile_list2_character <- map(shapefile_list2, convert_to_character)

# Combine the dataframes
combined_dataframe <- bind_rows(shapefile_list2_character)

################################################################################
#set data types

# Define the columns and their desired data types
numeric_columns <- c("long", "lat", "num_indepen", "num_pups_lar", "num_pups_sma", 
                     "point_x", "point_x")
factor_columns <- c("map_areas", "behavior", "habitat")

# Convert columns to their desired data types
combined_dataframe <- combined_dataframe %>%
  mutate(across(all_of(numeric_columns), as.numeric),
         across(all_of(factor_columns), as.factor),
         across(.cols = -c(all_of(numeric_columns), all_of(factor_columns)), as.character))


################################################################################
#clean up
census_dat_build1 <- combined_dataframe %>%
                      dplyr::select(-point_x, -point_y) %>%
                      #arrange
                      dplyr::select(long, lat, map_areas, date_time_orig = date_time, date_time = file_dt,
                                    observer1, observer2, num_indepen, num_pups_lar,
                                    num_pups_sma, everything()) %>%
                      #drop errors
                      filter(!(long %in% c(0.0000, "NaN"))) %>%
                      #make spatial with original geometry
                      st_as_sf(coords = c("long", "lat"), crs = 4326) 
    
#check
plot(census_dat_build1)                      


################################################################################
#save

output_path <- file.path("/Volumes/seaotterdb$/RESEARCH FORMS AND FILES/Monterey monthly census/QuarterlyCensus/qc_processed/merged_test_JGSv2.shp")

st_write(census_dat_build1, output_path, driver = "GeoJSON", append = TRUE)






