#Joshua G. Smith; jossmith@mbayaq.org

rm(list=ls())

######

#This script calculates the relative use of each scan area. 

######
#required packages
librarian::shelf(tidyverse, sf, raster, shiny, tmap, terra, tidyterra, ggnewscale)

#set directories 
basedir <- "/Volumes/seaotterdb$/kelp_recovery/data"

#read otter scan area
scan_area <- st_read(file.path(basedir, "gis_data/raw/otter_scan_area/TrackingMaps_POLY.shp")) 

# read census data from Step1
census_orig <- st_read(file.path("/Volumes/seaotterdb$/RESEARCH FORMS AND FILES/Monterey monthly census/QuarterlyCensus/qc_processed/merged_test_JGS.shp"))

################################################################################
# format scan area

scan_plot_area <- scan_area %>% mutate(Name = factor(Name)) %>%
  rename(site_name = Name)%>%
  mutate(Incipient = ifelse(site_name == "3200"|
                              site_name == "Carmel River Beach"|
                              site_name == "Coast Guard Pier" |
                              site_name == "El Torito" |
                              site_name == "Hopkins Marine Station East"|
                              site_name == "Lone Cypress" |
                              site_name =="Monastery" |
                              site_name == "Monterey Bay Inn" |
                              site_name == "Pescadero Point East" |
                              site_name == "Pescadero Point West"|
                              site_name == "Whaler's Cove","Yes","No"),
         incipient = factor(Incipient, levels = c("Yes","No"))) %>%
          dplyr::select(-Id, -Map_Number, -Incipient)


################################################################################
# process census data

census_build1 <- census_orig %>%
  #extract date
  mutate(
    year = lubridate::year(date_time),
    month = lubridate::month(date_time),
    day = lubridate::day(date_time),
    #set quarter
    quarter = case_when(
      month %in% c(1, 2, 3) ~ 1,
      month %in% c(4, 5, 6) ~ 2,
      month %in% c(7, 8, 9) ~ 3,
      month %in% c(10, 11, 12) ~ 4
    )
  ) %>% 
  #set geometry
  st_transform(crs = st_crs(scan_plot_area))%>%
  #join with scan area
  st_join(scan_plot_area)

#determine total counts by quarter and year for each behavior (foraging, resting, etc.)
census_build2 <- census_build1 %>%
  group_by(year, quarter, behavior)%>%
  dplyr::summarize(n_independ_total = sum(num_indepen),
                   n_pup_lar_total = sum(num_pups_lar),
                   n_pup_sm_total = sum(num_pups_sma)) 

#determine total counts by scan area and quarter for each behavior
census_build3 <- census_build1 %>%
  group_by(year, quarter, behavior, site_name)%>%
  dplyr::summarize(n_independ_site = sum(num_indepen),
                   n_pup_lar_site = sum(num_pups_lar),
                   n_pup_sm_site = sum(num_pups_sma)) %>%
  #join with total counts for the quarter
  st_join(census_build2) %>%
  #determine relative foraging
  mutate(n_rel_indepen = n_independ_site / n_independ_total,
         n_rel_sm = n_pup_sm_site / n_pup_sm_total,
         n_rel_lg = n_pup_lar_site / n_pup_lar_total
         )

#clean up
census_build4 <- census_build3 %>%
                  dplyr::select(year = year.x, quarter = quarter.x,
                                behavior = behavior.x, site_name, n_independ_site,
                                n_pup_lar_site, n_pup_sm_site, n_independ_total,
                                n_pup_lar_total, n_pup_sm_total, 
                                n_rel_indepen, n_rel_sm, n_rel_lg,
                                geometry
                                ) %>%
                  mutate(n_rel_sm = ifelse(n_rel_sm == "NaN",NA, n_rel_sm),
                         n_rel_lg = ifelse(n_rel_lg == "NaN",NA, n_rel_lg))


#join scan area
census_build5 <- scan_plot_area %>% st_join(census_build4) %>%
                  #clean up
                  dplyr::select(site_name = site_name.x,
                                incipient:behavior,
                                n_independ_site:geometry) %>%
                  #drop missing year and quarter
                  filter(!(is.na(year) | is.na(quarter)))












