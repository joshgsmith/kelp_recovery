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
envr_raw <- readRDS( file.path(basedir, "/data/environmental_data/cuti/processed/beuti_cuti_sst_monthly_amoms_by_PISCO_site.Rds"))
npp_raw <- readRDS( file.path(basedir, "/data/environmental_data/NPP/processed/npp_at_PISCO_site.Rds"))

#load environmental data from Natalie Low
nl_dat <- read.csv(file.path(basedir, "/data/environmental_data/PISCO_site_data/raw/Merged_Env_vars_all_sites_and_years.csv")) %>% janitor::clean_names()

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

envr_build1 <- envr_raw %>% #select sites in Carmel and Monterey Bay only
  left_join(site_table, by="site")%>%
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
  dplyr::select(!(c(latitude, longitude, ca_mpa_name_short, mpa_class, mpa_designation, baseline_region)))

################################################################################
#join substrate data with envr data

site_predict_build2 <- left_join(envr_build1, site_predict_build1, by="site") %>%
                      dplyr::select(year, month, baseline_region, site, latitude, longitude,
                                    ca_mpa_name_short, mpa_class, mpa_designation, everything()) %>%
                        #join with natalie low data
                        left_join(nl_dat, by=c("site","year")) %>%
                        #clean up
                      dplyr::select(!(c(id, campus, lat_wgs84, lon_wgs84, x_type, x_freq)))


#saveRDS(site_predict_build2, file.path(basedir, "data/environmental_data/predictors_at_pisco_sites_v2.Rds"))

################################################################################
#join with NPP

#npp_build1 <- npp_raw %>% dplyr::select(year, month, site, npp=chlorophyll) 

#site_predict_build3 <- site_predict_build2 %>% 
                        #rename to match missing values with their closest relative site
 #                       mutate(nearest_site = ifelse(site == "BLUEFISH_UC", "BLUEFISH_DC",
  #                                                                          ifelse(site == "BUTTERFLY_UC","BUTTERFLY_DC",
   #                                                                                ifelse(site == "PESCADERO_UC","PESCADERO_DC",site))))

#site_predict_build4 <- left_join(site_predict_build3, npp_build1, by=c("year","month","nearest_site"="site")) %>%
 #                       dplyr::select(!(nearest_site))



#export predictors

#saveRDS(site_predict_build4, file.path(basedir, "data/environmental_data/predictors_at_pisco_sites.Rds"))

################################################################################
#plot

library(ggplot2)


# create the plot
ggplot(npp_build1 %>% filter(npp <2), aes(x = year, y = npp, color = site
                                          )) +
  geom_point() +                    # add points
  stat_summary(fun.y = mean, geom = "line", size = 1.5) +   # add mean line
  scale_x_continuous(breaks = seq(min(npp_build1$year), max(npp_build1$year), by = 1)) +   # set breaks on x-axis
  labs(x = "Year", y = "NPP", color = "Site") +  # set axis and legend labels
  theme_bw()





ggplot(npp_build1 %>% filter(npp < 2), aes(x = as.Date(paste(year, month, "01", sep = "-")), y = npp, #color = site
                                           )) +
  geom_point(alpha=0.8, color="gray80", fill="gray80") +
  stat_summary(fun.y = mean, geom = "line", size = 1.5) +
  scale_x_date(date_breaks = "12 month", date_labels = "%Y", limits = c(as.Date("2003-01-01"), as.Date("2020-01-01"))) +
  labs(x = "Year", y = "NPP", color = "Site") +
  theme_bw()





