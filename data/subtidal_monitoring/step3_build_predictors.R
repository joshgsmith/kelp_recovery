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

#load environmental data
cuti_beuti_raw <- readRDS(file.path(basedir, "/data/environmental_data/CUTI/processed/1988_2022_cuti_beuti_daily_by_PISCO_site.Rds"))
sst_raw <- readRDS(file.path(basedir, "/data/environmental_data/MURSST/processed/2002_2022_mursst_monthly_by_PISCO_site.Rds"))

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
                    site == "BIRD_ROCK"|
                    site == "LINGCOD_DC"|
                    site == "LINGCOD_UC"|
                    site == "CARMEL_DC"|
                    site == "CARMEL_UC"))%>%
  distinct() #remove duplicates


################################################################################
#process seafloor data

vrm_build1 <- vrm_raw %>% dplyr::select(site = SITE_SIDE, vrm_mean = MEAN, 
                                        vrm_range = RANGE, vrm_sum = SUM)%>%
                #rename sites to match
  mutate(site = ifelse(site == "LONE_TREE_CEN","LONE_TREE", site),
         site = ifelse(site == "SIREN_CEN","SIREN", site))

slope_build1 <- slope_raw %>% dplyr::select(site = SITE_SIDE, slope_mean = MEAN, 
                                        slope_range = RANGE, slope_sum = SUM)%>%
  #rename sites to match
  mutate(site = ifelse(site == "LONE_TREE_CEN","LONE_TREE", site),
         site = ifelse(site == "SIREN_CEN","SIREN", site))

bat_build1 <- bat_raw %>% dplyr::select(site = SITE_SIDE, bat_mean = MEAN, 
                                            bat_range = RANGE) %>%
                        #make positive-definite
                        mutate(bat_mean = bat_mean *-1)%>%
  #rename sites to match
  mutate(site = ifelse(site == "LONE_TREE_CEN","LONE_TREE", site),
         site = ifelse(site == "SIREN_CEN","SIREN", site))


#join seafloor with site table

site_predict_build1 <- site_table %>% left_join(vrm_build1, by="site") %>%
                        left_join(slope_build1, by="site")%>%
                        left_join(bat_build1, by="site")


################################################################################
#process environmental data

cuti_beuti_build1 <- cuti_beuti_raw %>% dplyr::select(site, year, month, date, cuti, beuti) %>%
  #calculate annual means at each site
  group_by(site, year) %>%
  dplyr::summarize(cuti_avg = mean(cuti),
                   beuti_avg = mean(beuti))

sst_build1 <- sst_raw %>% 
  dplyr::mutate(year = lubridate::year(date),
                month = lubridate::month(date),
                day = lubridate::day(date))%>%
  dplyr::select(year, month, day, site, sst_c)%>%
  #calculate annial means at each site
  group_by(year, site)%>%
  dplyr::summarize(sst_c_mean = mean(sst_c))

################################################################################
#join sampling data with environmental

mod_dat_build1 <- left_join(swath_raw, cuti_beuti_build1, by = c("year", "site")) %>%
  dplyr::select(year, beuti_avg, cuti_avg, everything())

mod_dat_build2 <- left_join(mod_dat_build1, sst_build1, by = c("year", "site")) %>%
  dplyr::select(year, sst_c_mean, everything())

# create lagged variables
mod_dat_build3 <- mod_dat_build2 %>%
  mutate(lag_beuti_1 = lag(beuti_avg, 1),
         lag_cuti_1 = lag(cuti_avg, 1),
         lag_sst_1 = lag(sst_c_mean, 1),
         lag_beuti_2 = lag(beuti_avg, 2),
         lag_cuti_2 = lag(cuti_avg, 2),
         lag_sst_2 = lag(sst_c_mean, 2))



