#Communtiy change analyses
#Joshua G. Smith; jossmith@mbayaq.org

rm(list=ls())

librarian::shelf(tidyverse, here, ggplot2, mvabund)


################################################################################
#set directories and load data
basedir <- "/Volumes/seaotterdb$/kelp_recovery"
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

#load species attribute table
spp_attribute <- readxl::read_excel(file.path(basedir,"data/subtidal_monitoring/raw/spp_attribute_table.xlsx"), sheet = 2)

################################################################################
#process attribute table

spp_att_build1 <- spp_attribute %>% janitor::clean_names() %>% 
                    dplyr::select(genusspecies, common_name, primary_taxonomic, secondary_taxonomic,
                                  feeding_mode, primary_trophic) %>%
                    mutate(genusspecies =  str_to_sentence(gsub("_", " ", genusspecies))) %>%
                    #fix obsolete names
                    mutate(genusspecies_renamed = ifelse(genusspecies == "Megastrea gibberosa","Pomaulax gibberosus",
                                                         ifelse(genusspecies == "Parastichopus californicus","Apostichopus californicus",
                                                                ifelse(genusspecies == "Parastichopus parvimensis","Apostichopus parvimensis",
                                                                       ifelse(genusspecies == "Balanus nubilis","Balanus nubilus",
                                                                              ifelse(genusspecies == "Crassedoma giganteum","Crassadoma gigantea",
                                                                                     ifelse(genusspecies == "Loxorhynchus scyra spp","Loxorhynchus crispatus scyra acutifrons",
                                                                                            ifelse(genusspecies == "Pugettia spp","Pugettia foliata",
                                                                                                   ifelse(genusspecies == "Tethya aurantia","Tethya californiana",
                                                                                                          ifelse(genusspecies == "Coralline algae -crustose","Coralline algae crustose",
                                                                                                                 ifelse(genusspecies == "Red algae -encrusting","Red algae encrusting",
                                                                                                                        ifelse(genusspecies == "Red algae (leaf-like)","Red algae leaf like",
                                                                                                                               ifelse(genusspecies == "Serpulorbis squamigerus","Thylacodes squamigerus",
                                                                                                                                      ifelse(genusspecies == "Desmarestia spp","Desmarestia",
                                                                                                                                             ifelse(genusspecies == "Dictyotales spp","Dictyoneurum californicum reticulatum",
                                                                                                                                                    ifelse(genusspecies == "Sponge","Porifera",
                                                                                                                                                           ifelse(genusspecies == "Red algae (branching flat blade)","Red algae branching flat blade",
                                                                                                                                                                  ifelse(genusspecies == "Red algae (lacy branching)","Red algae lacy branching",
                                                                                                                                                                         ifelse(genusspecies == "Tunicate -colonial,compund,social","Tunicate colonial compund social",genusspecies)))))))))))))))))))
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
                  #dplyr::select(!(transition_site))%>%
                  #filter to significant species
                  left_join(swath_mvabund, by=c("transition_site","species"), relationship = "many-to-many")%>%
                  filter(!(is.na(group)))

upc_filtered <- upc_mod_dat %>%
  pivot_longer(12:50, names_to = "species", values_to = "counts")%>%
  #dplyr::select(!(transition_site))%>%
  #filter to significant species
  left_join(upc_mvabund, by=c("transition_site","species"), relationship = "many-to-many")%>%
  filter(!(is.na(group)))

fish_filtered <- fish_mod_dat %>%
  pivot_longer(12:67, names_to = "species", values_to = "counts")%>%
  #dplyr::select(!(transition_site))%>%
  #filter to significant species
  left_join(fish_mvabund, by=c("transition_site","species"), relationship = "many-to-many")%>%
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
  dplyr::summarize(mean_counts = mean(counts, na.rm=TRUE)+0.0001)%>% #ddd a small constant to avoid -Inf
  pivot_wider(names_from = outbreak_period,
              values_from = mean_counts) %>%
  mutate(perc_change = ((After - Before)/Before * 100))

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
  mutate(perc_change = ((After - Before)/Before * 100))

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
  mutate(perc_change = ((After - Before)/Before * 100))

# Calculate the average perc_change for each species
avg_perc_change <- fish_pc %>%
  group_by(species) %>%
  summarize(avg_change = mean(perc_change, na.rm = TRUE)) %>%
  arrange(desc(avg_change))

# Reorder the levels of the species factor based on avg_perc_change
fish_pc$species <- factor(fish_pc$species, levels = avg_perc_change$species)




##merge

plot_merge <- rbind(swath_pc, fish_pc, upc_pc) %>%
  mutate(species = str_to_sentence(gsub("_", " ", species)))%>%
  #join with species attributes
  left_join(spp_att_build1, by=c("species"="genusspecies_renamed")) %>%
  #drop sea stars
  filter(!(species == "Pycnopodia helianthoides" | species == "Orthasterias koehleri" |
             species == "Pisaster giganteus" | species == "Patiria miniata" | species == "Cirripedia")) %>%
  #rename common names
  mutate(common_name = ifelse(common_name == "Red Algae (Lacy Branching)","Red Algae (Lacy)",
                              ifelse(common_name == "Red Algae (Branching Flat Blade)", "Red algae (flat blade)",
                                     ifelse(common_name == "Decorator Crab, Moss Crab", "Decorator crab",
                                            ifelse(common_name == "Tunicate -Colonial,Compund,Social","Tunicate (colonial)",common_name)))),
         common_name = str_to_sentence(common_name),
         primary_trophic = str_to_sentence(primary_trophic))

avg_perc_change <- plot_merge %>%
  group_by(transition_site, species) %>%
  summarize(avg_change = mean(perc_change, na.rm = TRUE)) %>%
  arrange(desc(avg_change))



################################################################################
#Plot

my_theme <-  theme(axis.text=element_text(size=6),
                   axis.text.y = element_text(angle = 90, hjust = 0.5),
                   axis.title=element_text(size=8),
                   plot.tag=element_text(size=8, face="bold"),
                   plot.title =element_text(size=7, face="bold"),
                   # Gridlines 
                   panel.grid.major = element_blank(), 
                   panel.grid.minor = element_blank(),
                   panel.background = element_blank(), 
                   axis.line = element_line(colour = "black"),
                   # Legend
                   legend.key = element_blank(),
                   legend.background = element_rect(fill=alpha('blue', 0)),
                   legend.key.height = unit(1, "lines"), 
                   legend.text = element_text(size = 6),
                   legend.title = element_text(size = 7),
                   #legend.spacing.y = unit(0.75, "cm"),
                   #facets
                   strip.background = element_blank(),
                   strip.text = element_text(size = 6 ,face="bold"),
)

# Create ggplot plot

ggplot(plot_merge %>% group_by(transition_site) %>%
         mutate(species = tidytext::reorder_within(species, perc_change, transition_site),
                #make sentence case
                species = str_to_sentence(gsub("_", " ", species)),
                #remove trailing spaces
                species = str_remove(str_trim(species), "\\s+\\w+$")),
       aes(x = perc_change, y = reorder(species, perc_change, FUN = function(x) -max(x)))) +
  geom_point(aes(color = primary_trophic)) +
  geom_segment(aes(x = 0, xend = perc_change, yend = species, color = primary_trophic), linetype = "solid") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  #add genus species labels
  geom_text(
    aes(x = ifelse(plot_merge$perc_change >= 0, 0.1, -0.1), label = species, fontface = "italic"),
    hjust = ifelse(plot_merge$perc_change >= 0, 0, 1),
    color = "black",
    size = 3,
    position = position_nudge(y = -0.2)
  ) +
  #add common name labels
  geom_text(
    aes(x = ifelse(plot_merge$perc_change >= 0, 0.1, -0.1), label = common_name),
    hjust = ifelse(plot_merge$perc_change >= 0, 0, 1),
    color = "black",
    size = 3,
    position = position_nudge(y =0.2)
  ) +
  #add pre-post densities
  geom_text(
    aes(x = perc_change, label = paste0("(",round(Before,2),", ",round(After,2),")")),
    #hjust = ifelse(plot_merge$perc_change >= 0, 1, -1),
    color = "black",
    size = 3,
    position = position_nudge(y = 0.2)
  ) +
  scale_x_continuous(
    trans = ggallin::pseudolog10_trans,
    breaks = c(-20000, -10000, -1000, -100, -10, -1, 1, 10, 100, 1000, 10000),
    labels = c(-20000, -10000, -1000, -100, -10, -1, 1, 10, 100, 1000, 10000)
  ) +
  scale_color_brewer(palette = "Dark2")+
  facet_wrap(~ transition_site, ncol = 2, scales = "free_y") +
  xlab("Percentage Change") +
  ylab("Species") +
  theme_bw() + my_theme + theme(axis.text.y = element_blank())





resist_dat <- plot_merge %>% filter(transition_site == "no")

p1 <- ggplot(resist_dat,
       aes(x = perc_change, y = reorder(species, -perc_change))) +
  geom_point(aes(color = primary_trophic)) +
  geom_segment(aes(x = 0, xend = perc_change, yend = species, color = primary_trophic), linetype = "solid", size=1) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  #add genus species labels
  geom_text(
    aes(x = ifelse(resist_dat$perc_change >= 0, 0.1, -0.1), label = species, fontface = "italic"),
    hjust = ifelse(resist_dat$perc_change >= 0, 0, 1),
    color = "black",
    size = 2,
    position = position_nudge(y = -0.1)
  ) +
  #add common name labels
  geom_text(
    aes(x = ifelse(resist_dat$perc_change >= 0, 0.1, -0.1), label = common_name),
    hjust = ifelse(resist_dat$perc_change >= 0, 0, 1),
    color = "black",
    size = 3,
    position = position_nudge(y =0.1)
  ) +
  #add pre-post densities
  geom_text(
    aes(x = ifelse(resist_dat$perc_change >= 0, -7, 7), label = paste0("(", round(Before, 2), ", ", round(After, 2), ")")),
    color = "black",
    size = 3,
    # position = position_nudge(x = 0.15)
  ) +
  scale_x_continuous(
    trans = ggallin::pseudolog10_trans,
    breaks = c(-20000, -10000, -1000, -100, -10, -1, 1, 10, 100, 1000, 10000),
    labels = c(-20000, -10000, -1000, -100, -10, -1, 1, 10, 100, 1000, 10000)
  ) +
  scale_color_brewer(palette = "Dark2")+
  #facet_wrap(~ transition_site, ncol = 2, scales = "free_y") +
  xlab("") +
  ylab("") +
  labs(tag = "A", color = "Trophic role")+
  ggtitle("Kelp forest")+
  theme_bw() + my_theme + theme(axis.text.y = element_blank())
p1


transition_dat <- plot_merge %>% filter(transition_site == "yes") %>% filter(!(species == "Leptasterias hexactis" | species == "Cirripidia")) %>%
                    #drop UPC macro
                  filter(!(species == "Macrocystis pyrifera" & After <1))

p2 <- ggplot(transition_dat,
             aes(x = perc_change, y = reorder(species, -perc_change))) +
  geom_point(aes(color = primary_trophic)) +
  geom_segment(aes(x = 0, xend = perc_change, yend = species, color = primary_trophic), linetype = "solid", size=1) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  #add genus species labels
  geom_text(
    aes(x = ifelse(transition_dat$perc_change >= 0, 0.1, -0.1), label = species, fontface = "italic"),
    hjust = ifelse(transition_dat$perc_change >= 0, 0, 1),
    color = "black",
    size = 2,
    position = position_nudge(y = -0.2)
  ) +
  #add common name labels
  geom_text(
    aes(x = ifelse(transition_dat$perc_change >= 0, 0.1, -0.1), label = common_name),
    hjust = ifelse(transition_dat$perc_change >= 0, 0, 1),
    color = "black",
    size = 3,
    position = position_nudge(y =0.3)
  ) +
  # add pre-post densities
  geom_text(
    aes(x = ifelse(transition_dat$perc_change >= 0, -7, 7), label = paste0("(", round(Before, 2), ", ", round(After, 2), ")")),
    color = "black",
    size = 3,
   # position = position_nudge(x = 0.15)
  ) +
  scale_x_continuous(
    trans = ggallin::pseudolog10_trans,
    breaks = c(-10000, -1000, -100, -10, -1, 1, 10, 100, 1000, 10000),
    labels = c(-10000, -1000, -100, -10, -1, 1, 10, 100, 1000, 10000),
    limits = c(-10000, 20000)
  ) +
  scale_color_brewer(palette = "Dark2")+
  #facet_wrap(~ transition_site, ncol = 2, scales = "free_y") +
  xlab("") +
  ylab("") +
  labs(tag = "B", color = "Trophic role")+
  ggtitle("Sea urchin barren")+
  theme_bw() + my_theme + theme(axis.text.y = element_blank())
p2


combined_plot <- ggpubr::ggarrange(p1, p2, common.legend = TRUE, align = "h") 

combined_plot_annotated <- annotate_figure(combined_plot,
                bottom = text_grob("Percent change", 
                                   hjust = 4.5, x = 1, size = 10),
                left = text_grob("Species", rot = 90, size = 10)
)


ggsave(combined_plot_annotated, filename=file.path(figdir, "Fig4_mvabund.png"), bg = "white",
       width=8, height=10, units="in", dpi=600) 







