#Joshua G. Smith; jossmith@mbayaq.org

rm(list=ls())

######
#required packages
librarian::shelf(tidyverse, tidync, sf, raster, terra)

#set directories 
basedir <- "/Volumes/seaotterdb$/kelp_recovery/data"
zipdir <- "/Volumes/seaotterdb$/kelp_recovery/data/seafloor_data/raw/CSUMB_SFML/zip_files"
extract_path <- "/Volumes/seaotterdb$/kelp_recovery/data/seafloor_data/raw/CSUMB_SFML/extracted"
figdir <- here::here("analyses","figures")


################################################################################
#load data

# glance at files inside the archive
#MBsouth_central_zip <- unzip(file.path(zipdir, "MontereyBaySouthCentral/CMB_north_2m_habitat.zip"), list = TRUE)

#set the directory containing the .zip files
#zip_dir <- file.path(zipdir, "MontereyBaySouthCentral/zip")
zip_dir <- file.path(zipdir, "YankeePt")

# List all .zip files in the directory
zip_files <- list.files(zip_dir, pattern = ".zip$", full.names = TRUE)


# Loop through the .zip files and extract them to the 'extracted' folder
for (zip_file in zip_files) {
  # Extract the contents of the .zip file to a folder with the same name
  # as the .zip file in the "extracted" directory
  unzip(zip_file, exdir = file.path(paste0(extract_path,"/YankeePt"), gsub("\\.zip$", "", basename(zip_file))))
}



