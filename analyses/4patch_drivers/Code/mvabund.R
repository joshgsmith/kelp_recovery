#Communtiy change analyses
#Joshua G. Smith; jossmith@mbayaq.org

rm(list=ls())

librarian::shelf(tidyverse, here, ggplot2, mvabund)


################################################################################
#set directories and load data
basedir <- "/Volumes/seaotterdb$/kelp_recovery/"
figdir <- here::here("analyses","4patch_drivers","Figures")

#load raw dat
fish_raw <- read.csv(file.path(basedir, "data/subtidal_monitoring/processed/kelp_fish_counts_CC.csv")) %>%
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
                    site == "BIRD_ROCK")) 

upc_raw <- read.csv(file.path(basedir, "data/subtidal_monitoring/processed/kelp_upc_cov_CC.csv")) %>%
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
                    site == "BIRD_ROCK"))

swath_raw <- read.csv(file.path(basedir, "data/subtidal_monitoring/processed/kelp_swath_counts_CC.csv")) %>%
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
                    site == "BIRD_ROCK"))

################################################################################
#calculate transect means to reduce memory 

#drop species that were never encountered
fish_build1 <- fish_raw %>% dplyr::select(where(~ any(. != 0)))
upc_build1 <- upc_raw %>% dplyr::select(where(~ any(. != 0)))
swath_build1 <- swath_raw %>% dplyr::select(where(~ any(. != 0)))


################################################################################
#swath mvabund

swath_mod_dat <- swath_build1 %>% 
  mutate(outbreak_period = ifelse(year <2014, "Before","After")) %>%
  dplyr::select(outbreak_period, everything()) %>%
  dplyr::group_by(year, outbreak_period, MHW, baseline_region, latitude, longitude, site,
                  affiliated_mpa, mpa_class, mpa_designation) %>%
  dplyr::summarize(across(3:59, mean, na.rm = TRUE)) %>%
  #define transition sites
  mutate(transition_site = ifelse(site == "HOPKINS_UC" | site == "CANNERY_UC" |
                                    site == "SIREN" | site == "CANNERY_DC","no","yes"))%>%
  dplyr::select(transition_site, everything())

swath_transition <- swath_mod_dat %>% filter(transition_site == "yes")
swath_persist <- swath_mod_dat %>% filter(transition_site == "no")

#####run model for transition sites
#create multivariate object
swath_t_spp <- mvabund(swath_transition[, 12:68]) #exclude grouping vars
#fit the model
swath_t_model <- manyglm(swath_t_spp ~ swath_transition$outbreak_period)
#test for significance
swath_t_result <- anova.manyglm(swath_t_model, p.uni = "adjusted")
swath_t_out <- as.data.frame(swath_t_result[["uni.p"]])
#examine output
swath_t_sig <- swath_t_out %>%
  pivot_longer(cols=1:ncol(.), names_to="species")%>%
  drop_na()%>%
  filter(value <= 0.05) %>%
  mutate(group="swath",
         transition_site = "yes")


#####run model for persist sites
#create multivariate object
swath_p_spp <- mvabund(swath_persist[, 12:68]) #exclude grouping vars
#fit the model
swath_p_model <- manyglm(swath_p_spp ~ swath_persist$outbreak_period)
#test for significance
swath_p_result <- anova.manyglm(swath_p_model, p.uni = "adjusted")
swath_p_out <- as.data.frame(swath_p_result[["uni.p"]])
#examine output
swath_p_sig <- swath_p_out %>%
  pivot_longer(cols=1:ncol(.), names_to="species")%>%
  drop_na()%>%
  filter(value <= 0.05) %>%
  mutate(group="swath",
         transition_site = "no")

#merge

swath_mvabund <- rbind(swath_t_sig, swath_p_sig)


################################################################################
#upc mvabund


upc_mod_dat <- upc_build1 %>% 
  mutate(outbreak_period = ifelse(year <2014, "Before","After")) %>%
  dplyr::select(outbreak_period, everything()) %>%
  dplyr::group_by(year, outbreak_period, MHW, baseline_region, latitude, longitude, site,
                  affiliated_mpa, mpa_class, mpa_designation) %>%
  dplyr::summarize(across(3:41, mean, na.rm = TRUE)) %>%
  #define transition sites
  mutate(transition_site = ifelse(site == "HOPKINS_UC" | site == "CANNERY_UC" |
                                    site == "SIREN" | site == "CANNERY_DC","no","yes"))%>%
  dplyr::select(transition_site, everything())

upc_transition <- upc_mod_dat %>% filter(transition_site == "yes")
upc_persist <- upc_mod_dat %>% filter(transition_site == "no")

#####run model for transition sites
#create multivariate object
upc_t_spp <- mvabund(upc_transition[, 12:50]) #exclude grouping vars
#fit the model
upc_t_model <- manyglm(upc_t_spp ~ upc_transition$outbreak_period)
#test for significance
upc_t_result <- anova.manyglm(upc_t_model, p.uni = "adjusted")
upc_t_out <- as.data.frame(upc_t_result[["uni.p"]])
#examine output
upc_t_sig <- upc_t_out %>%
  pivot_longer(cols=1:ncol(.), names_to="species")%>%
  drop_na()%>%
  filter(value <= 0.05) %>%
  mutate(group="upc",
         transition_site = "yes")


#####run model for persist sites
#create multivariate object
upc_p_spp <- mvabund(upc_persist[, 12:50]) #exclude grouping vars
#fit the model
upc_p_model <- manyglm(upc_p_spp ~ upc_persist$outbreak_period)
#test for significance
upc_p_result <- anova.manyglm(upc_p_model, p.uni = "adjusted")
upc_p_out <- as.data.frame(upc_p_result[["uni.p"]])
#examine output
upc_p_sig <- upc_p_out %>%
  pivot_longer(cols=1:ncol(.), names_to="species")%>%
  drop_na()%>%
  filter(value <= 0.05) %>%
  mutate(group="upc",
         transition_site = "no")


#merge

upc_mvabund <- rbind(upc_t_sig, upc_p_sig)


################################################################################
#fish mvabund


fish_mod_dat <- fish_build1 %>% 
  mutate(outbreak_period = ifelse(year <2014, "Before","After")) %>%
  dplyr::select(outbreak_period, everything()) %>%
  dplyr::group_by(year, outbreak_period, MHW, baseline_region, latitude, longitude, site,
                  affiliated_mpa, mpa_class, mpa_designation) %>%
  dplyr::summarize(across(4:59, mean, na.rm = TRUE)) %>%
  #define transition sites
  mutate(transition_site = ifelse(site == "HOPKINS_UC" | site == "CANNERY_UC" |
                                    site == "SIREN" | site == "CANNERY_DC","no","yes"))%>%
  dplyr::select(transition_site, everything())

fish_transition <- fish_mod_dat %>% filter(transition_site == "yes")
fish_persist <- fish_mod_dat %>% filter(transition_site == "no")

#####run model for transition sites
#create multivariate object
fish_t_spp <- mvabund(fish_transition[, 12:67]) #exclude grouping vars
#fit the model
fish_t_model <- manyglm(fish_t_spp ~ fish_transition$outbreak_period)
#test for significance
fish_t_result <- anova.manyglm(fish_t_model, p.uni = "adjusted")
fish_t_out <- as.data.frame(fish_t_result[["uni.p"]])
#examine output
fish_t_sig <- fish_t_out %>%
  pivot_longer(cols=1:ncol(.), names_to="species")%>%
  drop_na()%>%
  filter(value <= 0.05) %>%
  mutate(group="fish",
         transition_site = "yes")


#####run model for persist sites
#create multivariate object
fish_p_spp <- mvabund(fish_persist[, 12:67]) #exclude grouping vars
#fit the model
fish_p_model <- manyglm(fish_p_spp ~ fish_persist$outbreak_period)
#test for significance
fish_p_result <- anova.manyglm(fish_p_model, p.uni = "adjusted")
fish_p_out <- as.data.frame(fish_p_result[["uni.p"]])
#examine output
fish_p_sig <- fish_p_out %>%
  pivot_longer(cols=1:ncol(.), names_to="species")%>%
  drop_na()%>%
  filter(value <= 0.05) %>%
  mutate(group="fish",
         transition_site = "no")


#merge

fish_mvabund <- rbind(fish_t_sig, fish_p_sig)



















