#Communtiy change analyses
#Joshua G. Smith; jossmith@mbayaq.org

rm(list=ls())

librarian::shelf(tidyverse, here, ggplot2, mvabund)


################################################################################
#set directories and load data
basedir <- "/Volumes/seaotterdb$/kelp_recovery"
figdir <- here::here("analyses","4patch_drivers","Figures")
tabdir <- here::here("analyses","4patch_drivers","Tables")

#load raw dat
swath_raw <- read.csv(file.path(basedir, "data/subtidal_monitoring/processed/kelp_swath_counts_CC.csv"))

upc_raw <- read.csv(file.path(basedir, "data/subtidal_monitoring/processed/kelp_upc_cov_CC.csv")) 

fish_raw <- read.csv(file.path(basedir, "data/subtidal_monitoring/processed/kelp_fish_counts_CC.csv"))

#load species attribute table
spp_attribute <- read.csv(file.path(tabdir,"TableS1_spp_table.csv")) %>% janitor::clean_names()

################################################################################
#calculate transect means to reduce memory 

#drop species that were never encountered
swath_build1 <- swath_raw %>% dplyr::select(where(~ any(. != 0)))
upc_build1 <- upc_raw %>% dplyr::select(where(~ any(. != 0)))
fish_build1 <- fish_raw %>% dplyr::select(where(~ any(. != 0)))


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
  mutate(outbreak_period = ifelse(year < 2014, "Before","After")) %>%
  dplyr::select(outbreak_period, everything()) %>%
  dplyr::group_by(year, outbreak_period, MHW, baseline_region, latitude, longitude, site,
                  affiliated_mpa, mpa_class, mpa_designation) %>%
  dplyr::summarize(across(4:57, mean, na.rm = TRUE)) %>%
  #define transition sites
  mutate(transition_site = ifelse(site == "HOPKINS_UC" | site == "CANNERY_UC" |
                                    site == "SIREN" | site == "CANNERY_DC","no","yes"))%>%
  dplyr::select(transition_site, everything())

fish_transition <- fish_mod_dat %>% filter(transition_site == "yes")
fish_persist <- fish_mod_dat %>% filter(transition_site == "no")

#####run model for transition sites
#create multivariate object
fish_t_spp <- mvabund(fish_transition[, 12:65]) #exclude grouping vars
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
fish_p_spp <- mvabund(fish_persist[, 12:65]) #exclude grouping vars
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
  pivot_longer(12:65, names_to = "species", values_to = "counts")%>%
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
  left_join(spp_attribute, by=c("species" = "genus_species")) %>%
  #drop sea stars
  filter(!(species == "Pycnopodia helianthoides" | species == "Orthasterias koehleri" |
             species == "Pisaster giganteus" | species == "Patiria miniata" | species == "Cirripedia" |
             species == "Mediaster aequalis")) %>%
  #rename common names
  mutate(common_name = str_to_sentence(common_name),
         trophic_ecology = str_to_sentence(trophic_ecology),
         #common_name = case_when(
          # common_name == "Red algae (lacy branching)" ~ "Red algae (lacy)",
          # common_name == "Red algae (branching flat blade)" ~ "Red algae (branching)",
          # common_name == "Decorator crab, moss crab" ~ "Decorator or moss crab",
          # common_name == "Tunicate colonial,compund,social" ~ "Colonial tunicate",
         #  TRUE ~common_name
        #   ),
         #fix trophic ecology
        # trophic_ecology = ifelse(common_name == "Red sea urchin","Herbivore",trophic_ecology),
        # species = case_when(
        #   species == "Red algae lacy branching" ~ "",
        #   species == "Red algae leaf like" ~ "",
           #species == "Coralline algae crustose" ~ "",
           #species == "Red algae encrusting" ~ "",
           #species == "Dictyoneurum californicum reticulatum" ~ "Dictyoneurum spp.",
           #species == "Tunicate colonial,compund,social" ~ "",
        #   TRUE ~ species
        # )
        # 
        ) %>%
  mutate(trophic_ecology = ifelse(trophic_ecology == "Autotroph","Primary producer",trophic_ecology))



avg_perc_change <- plot_merge %>%
  group_by(transition_site, species) %>%
  summarize(avg_change = mean(perc_change, na.rm = TRUE)) %>%
  arrange(desc(avg_change))



################################################################################
#Plot

my_theme <-  theme(axis.text=element_text(size=6, color = "black"),
                   axis.text.y = element_text(angle = 90, hjust = 0.5, color = "black"),
                   axis.title=element_text(size=8, color = "black"),
                   plot.tag=element_text(size=8, face="plain", color = "black"),
                   plot.title =element_text(size=7, face="bold", color = "black"),
                   # Gridlines 
                   panel.grid.major = element_blank(), 
                   panel.grid.minor = element_blank(),
                   panel.background = element_blank(), 
                   axis.line = element_line(colour = "black"),
                   # Legend
                   legend.key = element_blank(),
                   legend.background = element_rect(fill=alpha('blue', 0)),
                   legend.key.height = unit(1, "lines"), 
                   legend.text = element_text(size = 6, color = "black"),
                   legend.title = element_text(size = 7, color = "black"),
                   #legend.spacing.y = unit(0.75, "cm"),
                   #facets
                   strip.background = element_blank(),
                   strip.text = element_text(size = 6 ,face="bold", color = "black", hjust=0),
)

# Create ggplot plot

ggplot(plot_merge %>% group_by(transition_site) %>%
         mutate(species = tidytext::reorder_within(species, perc_change, transition_site),
                #make sentence case
                species = str_to_sentence(gsub("_", " ", species)),
                #remove trailing spaces
                species = str_remove(str_trim(species), "\\s+\\w+$")),
       aes(x = perc_change, y = reorder(species, perc_change, FUN = function(x) -max(x)))) +
  geom_point(aes(color = trophic_ecology)) +
  geom_segment(aes(x = 0, xend = perc_change, yend = species, color = trophic_ecology), linetype = "solid") +
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


# Define Dark2 color palette
dark2_palette <- c("#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E", "#E6AB02", "#A6761D", "#666666")


resist_dat <- plot_merge %>% filter(transition_site == "no") %>%
                mutate(trophic_ecology = ifelse(trophic_ecology == "Autotroph","Primary producer",trophic_ecology))
resist_dat$label <- with(resist_dat, ifelse(survey_method == "UPC", paste0("(", round(Before, 2), ", ", round(After, 2), ") *"), paste0("(", round(Before, 2), ", ", round(After, 2), ") \u2020")))

p1 <- ggplot(resist_dat,
             aes(x = perc_change, y = reorder(species, -perc_change))) +
  geom_point(aes(color = trophic_ecology)) +
  geom_segment(aes(x = 0, xend = perc_change, yend = species, color = trophic_ecology), linetype = "solid", size = 1) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  # add genus species labels
  geom_text(
    aes(x = ifelse(resist_dat$perc_change >= 0, 0.1, -0.1), label = species, fontface = "italic"),
    hjust = ifelse(resist_dat$perc_change >= 0, 0, 1),
    color = "black",
    size = 2,
    position = position_nudge(y = -0.1),
    data = resist_dat
  ) +
  # add common name labels
  geom_text(
    aes(x = ifelse(resist_dat$perc_change >= 0, 0.1, -0.1), label = common_name),
    hjust = ifelse(resist_dat$perc_change >= 0, 0, 1),
    color = "black",
    size = 3,
    position = position_nudge(y = 0.1),
    data = resist_dat
  ) +
  # add pre-post densities
  geom_text(
    aes(
      x = ifelse(resist_dat$perc_change >= 0, -8, 8),
      label = label
    ),
    color = "black",
    size = 3,
    data = resist_dat
  ) +
  scale_x_continuous(
    trans = ggallin::pseudolog10_trans,
    breaks = c(-20000, -10000, -1000, -100, -10, -1, 1, 10, 100, 1000, 10000),
    labels = c(-20000, -10000, -1000, -100, -10, -1, 1, 10, 100, 1000, 10000)
  ) +
  scale_color_manual(values = setNames(dark2_palette, unique(plot_merge$trophic_ecology))) +
  xlab("") +
  ylab("") +
  labs(tag = "A", color = "Trophic function") +
  ggtitle("Persistent") +
  theme_bw() +
  my_theme +
  theme(axis.text.y = element_blank())

p1


transition_dat <- plot_merge %>% filter(transition_site == "yes") %>% filter(!(species == "Leptasterias hexactis" | species == "Cirripidia")) %>%
                    #drop UPC macro
                   filter(!(survey_method == "UPC" & species == "Macrocystis pyrifera")) %>%
                  filter(!(survey_method == "Swath" & species == "Macrocystis pyrifera" & Before < 2)) %>%
  mutate(trophic_ecology = ifelse(trophic_ecology == "Autotroph","Primary producer",trophic_ecology))

# Create a new column for formatted labels
transition_dat$label <- with(transition_dat, ifelse(survey_method == "UPC", paste0("(", round(Before, 2), ", ", round(After, 2), ") *"), paste0("(", round(Before, 2), ", ", round(After, 2), ") \u2020" )))

p2 <- ggplot(transition_dat,
             aes(x = perc_change, y = reorder(species, -perc_change))) +
  geom_point(aes(color = trophic_ecology)) +
  geom_segment(aes(x = 0, xend = perc_change, yend = species, color = trophic_ecology), linetype = "solid", size = 1) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  # add genus species labels
  geom_text(
    aes(x = ifelse(transition_dat$perc_change >= 0, 0.1, -0.1), label = species, fontface = "italic"),
    hjust = ifelse(transition_dat$perc_change >= 0, 0, 1),
    color = "black",
    size = 2,
    position = position_nudge(y = -0.2)
  ) +
  # add common name labels
  geom_text(
    aes(x = ifelse(transition_dat$perc_change >= 0, 0.1, -0.1), label = common_name),
    hjust = ifelse(transition_dat$perc_change >= 0, 0, 1),
    color = "black",
    size = 3,
    position = position_nudge(y = 0.3)
  ) +
  # add pre-post densities
  geom_text(
    aes(
      x = ifelse(transition_dat$perc_change >= 0, -10, 10),
      label = label
    ),
    color = "black",
    size = 3
  ) +
  scale_x_continuous(
    trans = ggallin::pseudolog10_trans,
    breaks = c(-10000, -1000, -100, -10, -1, 1, 10, 100, 1000, 10000),
    labels = c(-10000, -1000, -100, -10, -1, 1, 10, 100, 1000, 10000),
    limits = c(-10000, 20000)
  ) +
  scale_color_manual(values = setNames(dark2_palette, unique(plot_merge$trophic_ecology))) +
  xlab("") +
  ylab("") +
  labs(tag = "B", color = "Trophic function") +
  ggtitle("Transitioned") +
  theme_bw() +
  my_theme +
  theme(axis.text.y = element_blank())

p2



combined_plot <- ggpubr::ggarrange(p1, p2, common.legend = TRUE, align = "h") 

combined_plot_annotated <- ggpubr::annotate_figure(combined_plot,
                bottom = ggpubr::text_grob("Percent change", 
                                   hjust = 4.6, x = 1, size = 10),
                left = ggpubr::text_grob("Species", rot = 90, size = 10, vjust=2)
)

combined_plot

ggsave(combined_plot_annotated, filename=file.path(figdir, "Fig4_mvabund2.png"), bg = "white",
       width=8.5, height=10, units="in", dpi=600) 







