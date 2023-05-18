#Communtiy change analyses
#Joshua G. Smith; jossmith@mbayaq.org

rm(list=ls())

librarian::shelf(tidyverse, here, ggplot2, mvabund)


################################################################################
#set directories and load data
basedir <- "/Volumes/seaotterdb$/kelp_recovery/"
figdir <- here::here("analyses","4patch_drivers","Figures")

#load raw dat
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
swath_build1 <- swath_raw %>% dplyr::select(where(~ any(. != 0)))
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
#swath mvabund


swath_mod_dat <- swath_build1 %>% 
  mutate(outbreak_period = ifelse(year <2014, "Before","After")) %>%
  dplyr::select(outbreak_period, everything()) %>%
  dplyr::group_by(year, outbreak_period, MHW, baseline_region, latitude, longitude, site,
                  affiliated_mpa, mpa_class, mpa_designation) %>%
  dplyr::summarize(across(4:59, mean, na.rm = TRUE)) %>%
  #define transition sites
  mutate(transition_site = ifelse(site == "HOPKINS_UC" | site == "CANNERY_UC" |
                                    site == "SIREN" | site == "CANNERY_DC","no","yes"))%>%
  dplyr::select(transition_site, everything())

swath_transition <- swath_mod_dat %>% filter(transition_site == "yes")
swath_persist <- swath_mod_dat %>% filter(transition_site == "no")

#####run model for transition sites
#create multivariate object
swath_t_spp <- mvabund(swath_transition[, 12:67]) #exclude grouping vars
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
swath_p_spp <- mvabund(swath_persist[, 12:67]) #exclude grouping vars
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
#filter data based on significant results

swath_filtered <- swath_mod_dat %>%
                    pivot_longer(12:68, names_to = "species", values_to = "counts")%>%
                  dplyr::select(!(transition_site))%>%
                  #filter to significant species
                  left_join(swath_mvabund, by="species", relationship = "many-to-many")%>%
                  filter(!(is.na(group)))

upc_filtered <- upc_mod_dat %>%
  pivot_longer(12:50, names_to = "species", values_to = "counts")%>%
  dplyr::select(!(transition_site))%>%
  #filter to significant species
  left_join(upc_mvabund, by="species", relationship = "many-to-many")%>%
  filter(!(is.na(group)))

swath_filtered <- swath_mod_dat %>%
  pivot_longer(12:67, names_to = "species", values_to = "counts")%>%
  dplyr::select(!(transition_site))%>%
  #filter to significant species
  left_join(swath_mvabund, by="species", relationship = "many-to-many")%>%
  filter(!(is.na(group)))



################################################################################
#plot using dumbbell approach

means_before <- swath_filtered %>%
  filter(outbreak_period == "Before") %>%
  group_by(species, transition_site) %>%
  summarize(mean_counts_before = mean(counts, na.rm = TRUE))

# Calculate the mean counts for "After" outbreak
means_after <- swath_filtered %>%
  filter(outbreak_period == "After") %>%
  group_by(species, transition_site) %>%
  summarize(mean_counts_after = mean(counts, na.rm = TRUE))

# Merge the mean counts for "Before" and "After"
means <- merge(means_before, means_after, by = c("species", "transition_site"))

# Create the dumbbell plot
ggplot(means) +
  geom_segment(aes(x = mean_counts_before, xend = mean_counts_after, y = species, yend = species, color = (mean_counts_after > mean_counts_before)), size = 1) +
  facet_wrap(~ transition_site, scales = "free_y") +
  scale_color_manual(values = c("blue", "red"), guide = FALSE) +
  labs(x = "Log(Counts)", y = "Species") +
  theme_minimal()



################################################################################
#prep data for plotting

# Calculate percent change for each species within each transition site
swath_pc <- swath_filtered %>%
  group_by(transition_site, outbreak_period, species)%>%
  dplyr::summarize(mean_counts = mean(counts, na.rm=TRUE))%>%
  pivot_wider(names_from = outbreak_period,
              values_from = mean_counts) %>%
  mutate(perc_change = ((After - Before)/After * 100))

# Calculate the average perc_change for each species
avg_perc_change <- swath_pc %>%
  group_by(species) %>%
  summarize(avg_change = mean(perc_change, na.rm = TRUE)) %>%
  arrange(desc(avg_change))

# Reorder the levels of the species factor based on avg_perc_change
swath_pc$species <- factor(swath_pc$species, levels = avg_perc_change$species)





# Calculate percent change for each species within each transition site
upc_pc <- upc_filtered %>%
  group_by(transition_site, outbreak_period, species)%>%
  dplyr::summarize(mean_counts = mean(counts, na.rm=TRUE))%>%
  pivot_wider(names_from = outbreak_period,
              values_from = mean_counts) %>%
  mutate(perc_change = ((After - Before)/After * 100))

# Calculate the average perc_change for each species
avg_perc_change <- upc_pc %>%
  group_by(species) %>%
  summarize(avg_change = mean(perc_change, na.rm = TRUE)) %>%
  arrange(desc(avg_change))

# Reorder the levels of the species factor based on avg_perc_change
upc_pc$species <- factor(upc_pc$species, levels = avg_perc_change$species)




# Calculate percent change for each species within each transition site
fish_pc <- fish_filtered %>%
  group_by(transition_site, outbreak_period, species)%>%
  dplyr::summarize(mean_counts = mean(counts, na.rm=TRUE))%>%
  pivot_wider(names_from = outbreak_period,
              values_from = mean_counts) %>%
  mutate(perc_change = ((After - Before)/After * 100))

# Calculate the average perc_change for each species
avg_perc_change <- fish_pc %>%
  group_by(species) %>%
  summarize(avg_change = mean(perc_change, na.rm = TRUE)) %>%
  arrange(desc(avg_change))

# Reorder the levels of the species factor based on avg_perc_change
fish_pc$species <- factor(fish_pc$species, levels = avg_perc_change$species)




##merge

plot_merge <- rbind(swath_pc, fish_pc, upc_pc)
avg_perc_change <- plot_merge %>%
  group_by(transition_site, species) %>%
  summarize(avg_change = mean(perc_change, na.rm = TRUE)) %>%
  arrange(desc(avg_change))


# Create ggplot plot
ggplot(plot_merge, aes(x = perc_change, y = reorder(species,-perc_change))) +
  geom_point() +
  geom_segment(aes(x = 0, xend = perc_change, yend = species), linetype = "solid") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  facet_wrap(~ transition_site, ncol = 2, scales = "free_y") +
  xlab("Percentage Change") +
  ylab("Species") +
  theme_bw()

















