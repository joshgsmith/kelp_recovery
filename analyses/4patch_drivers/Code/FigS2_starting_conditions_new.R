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

#calculate mean kelp density pre-2013

kelp_mean <- swath_sub1 %>% filter(species == "macrocystis_pyrifera") %>% group_by(year, site, outbreak_period)%>%
  summarise(kelp_mean = mean(counts, na.rm=TRUE),
            one_sd = sd(counts, na.rm=TRUE)) %>%
  #fix site names
  mutate(site = str_to_title(gsub("_", " ", site)),
         site = str_replace(site, "Dc", "DC"),
         site = str_replace(site, "Uc", "UC"))

################################################################################

# Perform t-test for each site
ttest_results <- kelp_mean %>%
  mutate(MHW = ifelse(year <2014, "before",ifelse(year > 2016, "after","during")))%>%
  filter(MHW == "before" | MHW == "after")%>%
  group_by(site) %>%
  summarise(p_value = t.test(kelp_mean ~ MHW)$p.value) %>%
  mutate(site = factor(site))

ttest_results$sig_col <- ifelse(ttest_results$p_value > 0.05, "forestgreen", "purple")

# Merge t-test results with kelp_mean data
kelp_mean <- left_join(kelp_mean, ttest_results, by = "site")


################################################################################
#plot

# Theme
my_theme <-  theme(axis.text=element_text(size=7),
                   axis.text.y = element_text(angle = 90, hjust = 0.5),
                   axis.title=element_text(size=9),
                   plot.tag=element_text(size=9),
                   plot.title =element_text(size=8, face="bold"),
                   # Gridlines 
                   panel.grid.major = element_blank(), 
                   panel.grid.minor = element_blank(),
                   panel.background = element_blank(), 
                   axis.line = element_line(colour = "black"),
                   # Legend
                   #legend.key = element_blank(),
                   legend.background = element_rect(fill=alpha('blue', 0)),
                   legend.key.height = unit(1, "lines"), 
                   legend.text = element_text(size = 10),
                   legend.title = element_text(size = 10),
                   #legend.spacing.y = unit(0.75, "cm"),
                   #facets
                   #strip.background = element_blank(),
                   strip.text = element_text(size = 6 ,face="bold"),
)


# Compute site-level mean and standard deviation of kelp_mean for the "before" outbreak period
before_kelp <- subset(kelp_mean, outbreak_period == "before")
site_stats <- before_kelp %>%
  group_by(site) %>%
  summarize(mean_kelp_before = mean(kelp_mean), sd_kelp_before = sd(kelp_mean)) %>%
  mutate(site = factor(site))


# Compute site-level mean and standard deviation of kelp_mean for the "after" outbreak period
after_kelp <- subset(kelp_mean, outbreak_period == "after")
site_stats_post <- after_kelp %>%
  group_by(site) %>%
  summarize(mean_kelp_post = mean(kelp_mean), sd_kelp_post = sd(kelp_mean)) %>%
  mutate(site = as.factor(site))

# Merge site-level mean and sd back into the kelp_mean dataset
kelp_mean1 <- kelp_mean %>% left_join(site_stats, by = "site")
kelp_mean <- kelp_mean1 %>% left_join(site_stats_post, by = "site")

# Create a new column indicating whether each year is within, above or below 1 sd of the mean at each site
kelp_mean$color <- ifelse(kelp_mean$kelp_mean > kelp_mean$mean_kelp_before + kelp_mean$sd_kelp_before, "Above 1 SD",
                          ifelse(kelp_mean$kelp_mean < kelp_mean$mean_kelp_before - kelp_mean$sd_kelp_before, "Below 1 SD", "Within 1 SD"))

# Subset the data to only include years exceeding +/- 1 sd of the mean at each site
subset_kelp <- kelp_mean[kelp_mean$color != "Within 1 SD",]


resist_sites <- kelp_mean %>% dplyr::filter(site =="Siren" | site == "Cannery DC" | site == "Hopkins UC" | site == "Cannery UC" 
                                            #| site == "Butterfly DC" #even though this site is non-sig, it went to 0. 
                                            ) %>% data.frame() %>% mutate(site = factor(site))
transition_sites <- kelp_mean %>% dplyr::filter(!(site =="Siren" | site == "Cannery DC" | site == "Hopkins UC" | site == "Cannery UC")) %>% data.frame()




#####################################

library(gridExtra)

# Sort resist_sites data frame by descending mean_kelp_before
resist_sites <- resist_sites[order(-resist_sites$mean_kelp_before), ]

# Determine the overall y-axis limits
y_limits <- with(resist_sites, range(kelp_mean - sd_kelp_before, kelp_mean + sd_kelp_before))

# Create a list to store individual plots
plot_list <- list()

# Loop through each site
for (site in resist_sites$site) {
  # Subset data for the current site
  site_data <- resist_sites[resist_sites$site == site, ]
  
  # Create the time series plot for the current site
  site_plot <- ggplot(data = site_data, aes(x = year, y = kelp_mean)) +
    geom_ribbon(aes(ymin = mean_kelp_before - sd_kelp_before, ymax = mean_kelp_before + sd_kelp_before), fill = "gray80", alpha = 0.5) +
    geom_line(show.legend = TRUE) +  # Include the legend
    #geom_point(data = subset_kelp %>% filter(site %in% site_data$site), aes(color = color), size = 2) +
    geom_segment(data = site_stats %>% filter(site %in% site_data$site),aes(x = 2007, xend = 2014, y = mean_kelp_before, yend = mean_kelp_before),linetype = "solid", color = "forestgreen") + # baseline
    geom_segment(data = site_stats %>% filter(site %in% site_data$site),aes(x = 2014, xend = 2020, y = mean_kelp_before, yend = mean_kelp_before),linetype = "dashed", color = "forestgreen") + # baseline
    geom_segment(data = site_stats_post %>% filter(site %in% site_data$site),aes(x = 2016, xend = 2020, y = mean_kelp_post, yend = mean_kelp_post),linetype = "solid", color = "purple") + # post
    labs(x = "", y = "", color = "", title = site) +
    scale_color_manual(values = c("Kelp density" = "black", "Baseline density" = "forestgreen", "Post-density" = "purple"),
                       labels = c("Kelp density", "Baseline density", "Post-density")) +
    scale_linetype_manual(values = c("solid", "dashed", "solid"),
                          labels = c("Kelp density", "Baseline density", "Post-density")) +
    scale_y_continuous(breaks = c(0, 100, 200, 300), limits = y_limits) +  # Set the y-axis limits
    theme_classic() +
    theme(plot.margin = unit(c(1, 0, 0, 0), "lines")) + my_theme +
    annotate(geom = "rect", xmin = 2014, xmax = 2016, ymin = -Inf, ymax = Inf, fill = "indianred1", alpha = 0.2) + my_theme +
    theme(legend.title = element_blank(),
          aspect.ratio = 1.1)  # Remove the legend title
  
  # Add the current site plot to the list
  plot_list[[site]] <- site_plot
}

# Combine all plots into a grid arrangement
grid_arrange <- do.call(grid.arrange, c(plot_list, ncol = 4))

persistent_plots <- annotate_figure(grid_arrange,
                              bottom = text_grob("Year", 
                                                 hjust = 13, vjust = -1, x = 1, size = 10),
                              left = text_grob("Kelp density \n(stipes per 60 m²)", 
                                                 hjust = 0.2, vjust = 0.2, x = 1, size = 10, rot=90)) + 
                    theme(panel.spacing = unit(0, "lines"),
                          plot.margin = margin(0, 0, 0, 0))

persistent_plots







# Sort transition_sites data frame by descending mean_kelp_before
transition_sites <- transition_sites[order(-transition_sites$mean_kelp_before), ]

# Determine the overall y-axis limits
y_limits <- with(transition_sites, range(kelp_mean - sd_kelp_before, kelp_mean + sd_kelp_before))

# Create a list to store individual plots
plot_list <- list()

# Loop through each site
for (site in transition_sites$site) {
  # Subset data for the current site
  site_data <- transition_sites[transition_sites$site == site, ]
  
  # Create the time series plot for the current site
  site_plot <- ggplot(data = site_data, aes(x = year, y = kelp_mean)) +
    geom_ribbon(aes(ymin = mean_kelp_before - sd_kelp_before, ymax = mean_kelp_before + sd_kelp_before), fill = "gray80", alpha = 0.5) +
    geom_line(show.legend = TRUE) +  # Include the legend
    #geom_point(data = subset_kelp %>% filter(site %in% site_data$site), aes(color = color), size = 2) +
    geom_segment(data = site_stats %>% filter(site %in% site_data$site),aes(x = 2007, xend = 2014, y = mean_kelp_before, yend = mean_kelp_before),linetype = "solid", color = "forestgreen") + # baseline
    geom_segment(data = site_stats %>% filter(site %in% site_data$site),aes(x = 2014, xend = 2020, y = mean_kelp_before, yend = mean_kelp_before),linetype = "dashed", color = "forestgreen") + # baseline
    geom_segment(data = site_stats_post %>% filter(site %in% site_data$site),aes(x = 2016, xend = 2020, y = mean_kelp_post, yend = mean_kelp_post),linetype = "solid", color = "purple") + # post
    labs(x = "", y = "", color = "", title = site) +
    scale_color_manual(values = c("Kelp density" = "black", "Baseline density" = "forestgreen", "Post-density" = "purple"),
                       labels = c("Kelp density", "Baseline density", "Post-density")) +
    scale_linetype_manual(values = c("solid", "dashed", "solid"),
                          labels = c("Kelp density", "Baseline density", "Post-density")) +
    scale_y_continuous(breaks = c(0, 100, 200, 300), limits = y_limits) +  # Set the y-axis limits
    theme_classic() +
    theme(plot.margin = unit(c(0, 0, 0, 0), "lines")) + my_theme +
    annotate(geom = "rect", xmin = 2014, xmax = 2016, ymin = -Inf, ymax = Inf, fill = "indianred1", alpha = 0.2) + my_theme +
    theme(legend.title = element_blank(),
          aspect.ratio = 1.1)  # Remove the legend title
  
  # Add the current site plot to the list
  plot_list[[site]] <- site_plot
}

# Combine all plots into a grid arrangement
grid_arrange <- do.call(grid.arrange, c(plot_list, ncol = 5))

transition_plots <- annotate_figure(grid_arrange,
                                    bottom = text_grob("Year", 
                                                       hjust = 13, vjust = -1, x = 1, size = 10),
                                    left = text_grob("Kelp density \n(stipes per 60 m²)", 
                                                     hjust = 0.2, vjust = 0.3, x = 1, size = 10, rot=90)) + 
  theme(panel.spacing = unit(0, "lines"),
        plot.margin = margin(0, 0, 0, 0))

transition_plots




g <- ggarrange(persistent_plots, transition_plots, ncol=1, heights = c(0.23, 0.77)
          ) 


ggsave(g, filename=file.path(figdir, "FigX_starting_conditions_new3.png"), 
       width=7.5, height=8, units="in", dpi=600, bg="white")





