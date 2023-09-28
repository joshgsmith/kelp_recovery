

#Joshua G. Smith; jossmith@mbayaq.org

rm(list=ls())


######
#required packages
librarian::shelf(tidyverse, sf)

#set directories 
basedir <- "/Volumes/seaotterdb$/kelp_recovery/data/census_data"
censusdir <- file.path(basedir,"annual_surveys/raw")
output <-file.path(basedir,"annual_surveys/processed")


################################################################################
#read raw Monterey counts 2019, 2021-23

# Get a list of all CSV files in the directory
csv_files <- list.files(path = file.path(censusdir,"usgs_mba_2019_2023_mtrypen"), 
                        pattern = "\\.csv$", full.names = TRUE)

# Initialize an empty list to store individual dataframes
df_list <- list()

# Loop through each CSV file, read it, and inspect the column names
for (csv_file in csv_files) {
  # Read the CSV file
  df <- read.csv(csv_file)
  
  # Print column names for inspection 
  cat("Column names for file", csv_file, ":\n")
  print(colnames(df))
  
  # Append the dataframe to the list
  df_list[[csv_file]] <- df
}

# NOTES AFTER INSPECTION
# 1. Files do not have the same column names. 
# 2. Some files are missing ATOS 

# Combine all dataframes into a single dataframe, filling NAs for missing columns
combined_df <- bind_rows(df_list, .id = "source_file")

# Reset row names 
rownames(combined_df) <- NULL


################################################################################
#clean


mpen_counts_build1 <- combined_df %>%
                        dplyr::select(year, longitude = lon, latitude = lat,
                                      ATOS, ind, smpup, lgpup, zone, habID)%>%
                        janitor::clean_names()


#Inspect

ggplot(mpen_counts_build1, aes(x = longitude, y = latitude)) +
  geom_point() +
  labs(title = "Scatter Plot of Points",
       x = "Longitude",
       y = "Latitude")


################################################################################
#save

saveRDS(mpen_counts_build1, file=file.path(output,"mpen_annual_counts_2019_23.Rds"))














