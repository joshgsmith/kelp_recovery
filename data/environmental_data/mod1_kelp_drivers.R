#Mixed effects model of patch drivers
#Joshua G. Smith
#May 9, 2023

rm(list=ls())
librarian::shelf(tidyverse, here, vegan, cowplot, ggpubr, lme4)

library(dplyr)

################################################################################
#set directories and load data
basedir <- "/Volumes/seaotterdb$/kelp_recovery/"
figdir <- here::here("analyses","4patch_drivers","Figures")

#load sampling data
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

#load environmental data

cuti_beuti_raw <- readRDS(file.path(basedir, "/data/environmental_data/CUTI/processed/1988_2022_cuti_beuti_daily_by_PISCO_site.Rds"))
sst_raw <- readRDS(file.path(basedir, "/data/environmental_data/MURSST/processed/2002_2022_mursst_monthly_by_PISCO_site.Rds"))

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

################################################################################
#check for collinearity and normalize where needed

# Select the predictor variables to check for collinearity
predictors <- c("beuti_avg", "cuti_avg", "sst_c_mean")

# Calculate the correlation matrix between the predictor variables
cor_matrix <- cor(mod_dat_build3[predictors])

# Print the correlation matrix
print(cor_matrix)

#normalize predictors to account for collinearity
mod_dat_build3$sst_norm <- scale(mod_dat_build3$sst_c_mean)

################################################################################
#build model

# Fit mixed model with random intercepts and slopes for site and year
mod1 <- lmer(macrocystis_pyrifera ~ sst_c_mean + lag_sst_1 + #lag_sst_2+#lag_beuti_1 + lag_cuti_1 +#set fixed effects
               (1 | year:site),  #set random effects grouping
             data = mod_dat_build3)

# Print model summary
summary(mod1)

# Add residuals to the data frame
# Create partial residual plots for each predictor
car::crPlots(mod1)

# Create facet plot of residuals by site
ggplot(data = mod_dat_build3, aes(x = macrocystis_pyrifera, y = resid, color = site)) +
  geom_point() +
  facet_wrap(~site, ncol = 4) +
  labs(x = "Observed Macrocyctis Pyrifera", y = "Residuals") +
  theme_minimal()

















