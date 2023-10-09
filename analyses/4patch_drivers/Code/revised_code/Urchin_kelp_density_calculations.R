#Community data processing
#Joshua G. Smith
#May 1, 2023

rm(list=ls())
librarian::shelf(tidyverse, here, vegan)


################################################################################
#set directories and load data
basedir <- "/Volumes/seaotterdb$/kelp_recovery/"
figdir <- here::here("analyses","4patch_drivers","Figures")

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
#select and reshape data

swath_sub1 <- swath_raw %>% dplyr::select(1:11, 'macrocystis_pyrifera', 'strongylocentrotus_purpuratus') %>%
  pivot_longer(cols = c('macrocystis_pyrifera', 'strongylocentrotus_purpuratus'),
               names_to = "species", values_to="counts") %>%
  mutate(outbreak_period = ifelse(year <2014, "before","after"))%>%
  dplyr::select(year, outbreak_period, everything())%>%
  data.frame() 

#calculate mean transect density before vs. after
mean_density <- swath_sub1 %>% 
  group_by(outbreak_period, species)%>%
  summarise(mean_density_60m = mean(counts, na.rm=TRUE),
            one_sd_60m = sd(counts, na.rm=TRUE),
            mean_density_m2 = mean_density_60m/60,
            one_sd_m = one_sd_60m/60
            ) 

percent_loss <- mean_density %>%
  filter(outbreak_period %in% c("before", "after") & species == "macrocystis_pyrifera") %>%
  group_by(outbreak_period) %>%
  summarise(mean_density = mean(mean_density_60m)) %>%
  spread(outbreak_period, mean_density) %>%
  mutate(percent_loss = ((before - after) / before) * 100) %>%
  pull(percent_loss)

#calculate mean transect density before vs. 2020
mean_density <- swath_sub1 %>% 
  mutate(final_year = ifelse(year < 2014, "before",
                             ifelse(year == 2020,"after",NA)))%>%
  dplyr::filter(!(is.na(final_year)))%>%
  group_by(outbreak_period, species)%>%
  summarise(mean_density_60m = mean(counts, na.rm=TRUE),
            one_sd_60m = sd(counts, na.rm=TRUE),
            mean_density_m2 = mean_density_60m/60,
            one_sd_m = one_sd_60m/60
  ) 


percent_loss_2020 <- mean_density %>%
  filter(outbreak_period %in% c("before", "after") & species == "macrocystis_pyrifera") %>%
  group_by(outbreak_period) %>%
  summarise(mean_density = mean(mean_density_60m)) %>%
  spread(outbreak_period, mean_density) %>%
  mutate(percent_loss = ((before - after) / before) * 100) %>%
  pull(percent_loss)


################################################################################
#test using ANOVA

mod <- aov(counts~ species * outbreak_period * site, data = swath_sub1)
summary(mod)

TukeyHSD(mod)


anova_table <- summary(mod)
# Extract the numerator, denominator, and F-value
numerator <- anova_table[[1]]$"F value"[1]
denominator <- anova_table[[1]]$"Df"[2]
F_value <- anova_table[[1]]$"F value"[2]





