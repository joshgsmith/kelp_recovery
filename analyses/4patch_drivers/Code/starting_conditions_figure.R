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
                   legend.key = element_blank(),
                   legend.background = element_rect(fill=alpha('blue', 0)),
                   legend.key.height = unit(1, "lines"), 
                   legend.text = element_text(size = 10),
                   legend.title = element_text(size = 10),
                   #legend.spacing.y = unit(0.75, "cm"),
                   #facets
                   #strip.background = element_blank(),
                   strip.text = element_text(size = 6 ,face="bold"),
)



#######THIS WORKS GOOD

library(dplyr)
library(ggplot2)

# Compute site-level mean and standard deviation of kelp_mean for the "before" outbreak period
before_kelp <- subset(kelp_mean, outbreak_period == "before")
site_stats <- before_kelp %>%
  group_by(site) %>%
  summarize(mean_kelp = mean(kelp_mean), sd_kelp = sd(kelp_mean))

# Merge site-level mean and sd back into the kelp_mean dataset
kelp_mean <- kelp_mean %>% left_join(site_stats, by = "site")

# Create a new column indicating whether each year is within, above or below 1 sd of the mean at each site
kelp_mean$color <- ifelse(kelp_mean$kelp_mean > kelp_mean$mean_kelp + kelp_mean$sd_kelp, "Above 1 SD",
                          ifelse(kelp_mean$kelp_mean < kelp_mean$mean_kelp - kelp_mean$sd_kelp, "Below 1 SD", "Within 1 SD"))

# Subset the data to only include years exceeding +/- 1 sd of the mean at each site
subset_kelp <- kelp_mean[kelp_mean$color != "Within 1 SD",]


resist_sites <- kelp_mean %>% dplyr::filter(site == "Hopkins UC" | site == "Cannery UC" | site == "Macabee DC" | site =="Siren" | site == "Cannery DC") %>% data.frame()
transition_sites <- kelp_mean %>% dplyr::filter(!(site == "Hopkins UC" | site == "Cannery UC" | site == "Macabee DC" | site =="Siren" | site == "Cannery DC")) %>% data.frame()


# Create the time series plot
overall <- ggplot(data = kelp_mean, aes(x = year, y = kelp_mean)) +
  geom_line() +
  geom_ribbon(data = kelp_mean , aes(ymin = mean_kelp - sd_kelp, ymax = mean_kelp + sd_kelp), fill = "gray80", alpha = 0.5) +
  geom_point(data = subset_kelp, aes(color = color), size = 3) +
  geom_hline(data = site_stats , aes(yintercept = mean_kelp), linetype = "solid") + # add horizontal line
  labs(x = "Year", y = "Kelp Mean", color = "") +
  scale_color_manual(values = c("Within 1 SD" = "black", "Above 1 SD" = "forestgreen", "Below 1 SD" = "purple")) +
  facet_wrap(~reorder(site, -mean_kelp), scales = "fixed", ncol=5) +
  theme_classic()

overall



####split

# Create the time series plot for resist_sites
panel_A <- ggplot(data = resist_sites, aes(x = year, y = kelp_mean)) +
  geom_line() +
  geom_ribbon(aes(ymin = mean_kelp - sd_kelp, ymax = mean_kelp + sd_kelp), fill = "gray80", alpha = 0.5) +
  geom_point(data = subset_kelp %>% filter(site %in% resist_sites$site), aes(color = color), size = 2) +
  geom_hline(data = site_stats %>% filter(site %in% resist_sites$site), aes(yintercept = mean_kelp), linetype = "solid") + # add horizontal line
  labs(x = "Year", y = "Kelp density \n(stipes per 60 m²)", color = "") +
  scale_color_manual(values = c("Within 1 SD" = "black", "Above 1 SD" = "forestgreen", "Below 1 SD" = "purple")) +
  facet_wrap(~reorder(site, -mean_kelp), scales = "fixed", nrow=1, ncol=5) +
  scale_y_continuous(breaks = c(0, 100,200,300))+
  labs(tag = "A") +
  theme_classic() +
  theme(plot.margin = unit(c(0, 1, 0.5, 0), "lines"),
        aspect.ratio = 1.1) + my_theme

# Create the time series plot for transition_sites
panel_B <- ggplot(data = transition_sites, aes(x = year, y = kelp_mean)) +
  geom_line() +
  geom_ribbon(aes(ymin = mean_kelp - sd_kelp, ymax = mean_kelp + sd_kelp), fill = "gray80", alpha = 0.5) +
  geom_point(data = subset_kelp %>% filter(site %in% transition_sites$site), aes(color = color), size = 2) +
  geom_hline(data = site_stats %>% filter(site %in% transition_sites$site), aes(yintercept = mean_kelp), linetype = "solid") + # add horizontal line
  labs(x = "Year", y = "Kelp density \n(stipes per 60 m²)", color = "") +
  scale_color_manual(values = c("Within 1 SD" = "black", "Above 1 SD" = "forestgreen", "Below 1 SD" = "purple")) +
  facet_wrap(~reorder(site, -mean_kelp), scales = "fixed", ncol=5) +
  scale_y_continuous(breaks = c(0, 100,200,300), limits = c(0,300))+
  labs(tag = "B") +
  theme_classic() +
  theme(plot.margin = unit(c(0, 1, 0.5, 0), "lines"),
        aspect.ratio = 1.1) + my_theme

# Combine the two plots into one figure with two panels
g <- ggpubr::ggarrange(panel_A, panel_B, ncol = 1, common.legend = TRUE, legend = "top",
                  align = "v", heights = c(1, 3.25)
                  ) # + theme(panel.spacing = unit(0.01, "lines"))


g

ggsave(g, filename=file.path(figdir, "FigX_starting_conditions.png"), 
       width=7, height=8, units="in", dpi=600, bg="white")





