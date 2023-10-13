

#Joshua G. Smith; jossmith@mbayaq.org

rm(list=ls())


######
#required packages

librarian::shelf(tidyverse, sf, raster, terra, janitor)

basedir <- "/Volumes/seaotterdb$/kelp_recovery/data"
figdir <- here::here("analyses","5community_regulation","figures")


#read sofa output
mass_class <- readxl::read_excel(file.path(basedir,"/sofa_data/raw/MCResults_Periods_10-11-23_12h.xlsx"), sheet = 4, skip=1) %>% clean_names()
dietcomp_class <- readxl::read_excel(file.path(basedir,"/sofa_data/raw/MCResults_Periods_10-11-23_12h.xlsx"), sheet = 5, skip = 1) %>% clean_names()
kcal_class <- readxl::read_excel(file.path(basedir,"/sofa_data/raw/MCResults_Periods_10-11-23_12h.xlsx"), sheet = 6, skip=1) %>% clean_names()



################################################################################
#Step 1 - process sofa

ncol(mass_class)
ncol(dietcomp_class)
ncol(kcal_class)

# Add a column to each dataset indicating the source
mass_class1 <- mass_class %>% mutate(source = "mass_gain")
dietcomp_class1 <- dietcomp_class %>% mutate(source = "dietcomp")
kcal_class1 <- kcal_class %>% mutate(source = "kcal")

# Combine the datasets into one
sofa_build1 <- bind_rows(mass_class1, dietcomp_class1, kcal_class1)

################################################################################
#Step 2 - process intake
kcal_sum <- kcal_class %>%
  filter(period <=2019)%>%
  rename(year = period) %>%
  mutate(period = ifelse(year < 2013, "Pre-SSW", "Post-SSW")) %>%
  group_by(period) %>%
  summarize(across(.cols = 2:25, .fns = mean, na.rm = TRUE)) %>%
  pivot_longer(cols = 2:25, names_to = "prey", values_to = "kcal") %>%
  filter(!prey %in% c("lobster", "fish", "fishegg", "small_kelp_invert")) %>%
  mutate(predator = "Sea otter",
         period = factor(period, levels = c("Pre-SSW","Post-SSW")),
         prey_kcal = paste(prey,round(kcal,2))) %>%
  filter(!(kcal<0.01))

hist(kcal_sum$kcal)

# Calculate percent change and format prey
kcal_sum_perc <- kcal_sum %>%
  dplyr::select(period, prey, kcal) %>%
  pivot_wider(names_from = period, values_from = kcal) %>%
  filter(!(prey %in% c('sand_crab'))) %>%
  mutate(
    prey = str_replace_all(prey, "_", " "),  # Remove underscores
    prey = tolower(prey),                    # Convert to lowercase
    prey = str_to_sentence(prey),             # Convert to sentence case
    perc_change = ((`Post-SSW` - `Pre-SSW`) / `Pre-SSW`) * 100,
    label = paste0("(", round(`Pre-SSW`, 2), ", ", round(`Post-SSW`, 2), ")"),
    #assign functional group
    functional_group = case_when(
      prey == "Mussel" ~ "Planktivore",
      prey == "Urchin" ~ "Grazer",
      prey == "Snail" ~ "Grazer",
      prey == "Abalone" ~ "Grazer",
      prey == "Clam" ~ "Planktivore",
      prey == "Octopus" ~ "Macroinvertivore",
      prey == "Squid" ~ "Macroinvertivore",
      prey == "Kelp crab" ~ "Grazer",
      prey == "Worm" ~ "Detritivore",
      prey == "Cancer crab" ~ "Macroinvertivore",
      prey == "Sand dollar" ~ "Planktivore",
      prey == "Crab other" ~ "Macroinvertivore",
      prey == "Chiton" ~ "Grazer",
      prey == "Star" ~ "Macroinvertivore"
    )
  )



################################################################################
#Step 4 - plot


# Theme
my_theme <-  theme(axis.text=element_text(size=8,color = "black"),
                   axis.title=element_text(size=8,color = "black"),
                   plot.tag=element_text(size=7,color = "black"),
                   plot.title=element_text(size=10,color = "black", face = "bold"),
                   # Gridlines
                   panel.grid.major = element_blank(), 
                   panel.grid.minor = element_blank(),
                   panel.background = element_blank(), 
                   axis.line = element_line(colour = "black"),
                   # Legend
                   legend.key.size = unit(0.3, "cm"), 
                   #legend.key = element_rect(fill = "white"), # Set it to transparent
                   legend.spacing.y = unit(0.1, "cm"),  
                   legend.text=element_text(size=8,color = "black"),
                   legend.title=element_blank(),
                   #legend.key.height = unit(0.1, "cm"),
                   #legend.background = element_rect(fill=alpha('blue', 0)),
                   #facets
                   strip.text = element_text(size=8, face = "bold",color = "black", hjust=0),
                   strip.background = element_blank())


# Reshape the data from wide to long format
stacked_data <- sofa_build1 %>%
  dplyr::select(period, source, 2:22) %>%
  pivot_longer(cols = -c(period, source), names_to = "Prey", values_to = "Value") 

# Filter the prey types that reach or exceed 3 for mass_gain
prey_to_include <- stacked_data %>%
  filter(source == "mass_gain", Value >= 3) %>%
  distinct(Prey)

# Filter the data to include the selected prey types
filtered_data <- stacked_data %>%
  filter(Prey %in% c("urchin","mussel"))


a <- ggplot(filtered_data %>% filter(period < 2020),
       aes(x = period, y = Value, color = Prey)) +
  #geom_line() +
  geom_point(alpha = 0.4) +
  #add SSW reference
  geom_vline(xintercept = 2013, linetype = "dashed", color = "black")+
  #geom_smooth(method = "loess", se = TRUE, aes(fill = Prey)) + 
  stat_smooth(geom = "line", size = 1, span = 0.7) +
  stat_smooth(method = "loess", geom = "ribbon", alpha = 0.2, aes(fill = Prey), color = NA) +
  facet_wrap(~ source, ncol = 1, scales = "free_y", 
             labeller = labeller(source = 
                                   c(dietcomp = "Proportion diet",
                                     kcal = "Energy gain (Kcal / min)",
                                     mass_gain = "Mass gain (g/min)"))) +
  labs(
    title = "",
    x = "Year",
    y = "Prey value",
    tag = "A"
  ) +
  scale_x_continuous(breaks = seq(2007, 2020, by = 2)) +  
  theme_bw() + my_theme+
  theme(legend.position = "top")+
  scale_color_brewer(palette = "Dark2")+
  scale_fill_brewer(palette = "Dark2") 
  #theme(plot.tag.position = c(0,0.88))
a




# Create a dumbbell plot
b <- ggplot(kcal_sum_perc, aes(x = perc_change, y = reorder(prey, perc_change)
)) +
  geom_point(aes(color = functional_group, fill = functional_group)) +
  geom_segment(aes(x = 0, xend = perc_change, yend = prey,color = functional_group), linetype = "solid", size = 1) +
  geom_text(aes( x = ifelse(kcal_sum_perc$perc_change >= 0, -50, 50),
                 label = label), size=2.5) +  
  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(x = "Percent change", y = "", tag = "B") +
  scale_x_continuous(limits=c(-100,350))+
  theme_bw() + my_theme + theme(legend.position = "top")+
  guides(colour = guide_legend(nrow = 2))
b


g <- gridExtra::grid.arrange(a,b, ncol=2)

g

ggsave(g, filename = file.path(figdir, "Fig3_foraging_timeseriesv2.png"), 
      width = 7, height = 6, units = "in", dpi = 600, bg = "white")


