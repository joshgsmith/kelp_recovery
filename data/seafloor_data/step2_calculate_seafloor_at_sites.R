##Joshua G. Smith
#May 10, 2023

rm(list=ls())
librarian::shelf(tidyverse, here, readxl)


#rasters were obtained from the SeaFloor Mapping Lab at CSU Monterey Bay
#and processed in ArcGIS. The following script links the processed tables 
#to the subtidal data

################################################################################
#set directories and load data
basedir <- "/Volumes/seaotterdb$/kelp_recovery"
figdir <- here::here("analyses","4patch_drivers","Figures")

#load monitoring site table
site_table <- read.csv(file.path(basedir, "data/subtidal_monitoring/raw/MLPA_kelpforest_site_table.4.csv")) %>%
  janitor::clean_names() %>%
  dplyr::select(site, latitude, longitude, ca_mpa_name_short, mpa_class=site_designation, mpa_designation=site_status, 
                baseline_region)%>%
  distinct()%>% #remove duplicates
  #filter to central coast region
  filter(baseline_region == "CENTRAL")

#load seafloor data

vrm_raw <- read_excel(file.path(basedir, "data/seafloor_data/processed/vrm_pisco_sites.xls"))
bathy_raw <- read_excel(file.path(basedir, "data/seafloor_data/processed/bathy_pisco_sites.xls"))
slope_raw <- read_excel(file.path(basedir, "data/seafloor_data/processed/slope_pisco_sites.xls"))

################################################################################
#process and merge seafloor metrics

vrm_build1 <- vrm_raw %>% dplyr::select(!(c(OID, ZONE_CODE, COUNT, AREA))) 
names(vrm_build1)[2:7] <- paste0("vrm_", names(vrm_build1)[2:7])
                  
bathy_build1 <- bathy_raw %>% dplyr::select(!(c(Rowid, `ZONE-CODE`, COUNT, AREA))) 
names(bathy_build1)[2:7] <- paste0("bathy_", names(bathy_build1)[2:7])

slope_build1 <- slope_raw %>% dplyr::select(!(c(Rowid, `ZONE-CODE`, COUNT, AREA))) 
names(slope_build1)[2:7] <- paste0("slope_", names(slope_build1)[2:7])

#join
joined_df <- inner_join(vrm_build1, bathy_build1, by = "SITE_SIDE") %>%
  inner_join(slope_build1, by = "SITE_SIDE")

################################################################################
#process and merge seafloor metrics

site_seafloor <- left_join(site_table, joined_df, by = c("site"= "SITE_SIDE"))%>%
                  filter(!is.na(vrm_MEAN))


write.csv(site_seafloor, file.path(basedir, "data/seafloor_data/processed/site_seafloor_data.csv"), row.names = FALSE)







