##Joshua G. Smith
#May 15, 2023

rm(list=ls())
librarian::shelf(tidyverse)


# Set Directories
basedir <-"/Volumes/seaotterdb$/kelp_recovery/"
figdir <- here::here("analyses","4patch_drivers","Figures")

#load seafloor data
vrm_raw <- readxl::read_excel(file.path(basedir,"data/seafloor_data/processed/vrm_pisco_sites.xls"))
slope_raw  <- readxl::read_excel(file.path(basedir,"data/seafloor_data/processed/slope_pisco_sites.xls"))
bat_raw <- readxl::read_excel(file.path(basedir,"data/seafloor_data/processed/bathy_pisco_sites.xls"))

################################################################################
#load sites

#load monitoring site table
site_table <- read.csv(file.path(basedir, "data/subtidal_monitoring/raw/MLPA_kelpforest_site_table.4.csv")) %>%
  janitor::clean_names() %>%
  dplyr::select(site, latitude, longitude, ca_mpa_name_short, mpa_class=site_designation, mpa_designation=site_status, 
                baseline_region)%>%
  #select central coast as focal region
  #select sites in Carmel and Monterey Bay only
  dplyr::filter(latitude >= 36.46575 & latitude <= 36.64045) %>%
  #drop sites with insufficient data
  dplyr::filter(!(site == "ASILOMAR_DC" |
                    site == "ASILOMAR_UC" |
                    site == "CHINA_ROCK" |
                    site == "CYPRESS_PT_DC" |
                    site == "CYPRESS_PT_UC" |
                    site == "PINNACLES_IN" |
                    site == "PINNACLES_OUT" |
                    site == "PT_JOE" |
                    site == "SPANISH_BAY_DC" |
                    site == "SPANISH_BAY_UC" |
                    site == "BIRD_ROCK"))%>%
  distinct() #remove duplicates


################################################################################
#process VRM






