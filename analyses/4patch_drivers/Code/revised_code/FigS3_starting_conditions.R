#Community data processing
#Joshua G. Smith
#May 1, 2023

rm(list=ls())
librarian::shelf(tidyverse, here, vegan)


################################################################################
#set directories and load data
basedir <- "/Volumes/seaotterdb$/kelp_recovery/"
figdir <- here::here("analyses","4patch_drivers","Figures")

swath_raw <- read.csv(file.path(basedir, "data/subtidal_monitoring/processed/kelp_swath_counts_CC.csv")) 

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


resist_sites <- kelp_mean %>% dplyr::filter(site =="Siren" | site == "Cannery DC" | site == "Hopkins UC" | site == "Cannery UC" 
                                            #| site == "Butterfly DC" #even though this site is non-sig, it went to 0. 
                                            ) %>% data.frame()
transition_sites <- kelp_mean %>% dplyr::filter(!(site =="Siren" | site == "Cannery DC" | site == "Hopkins UC" | site == "Cannery UC")) %>% data.frame()


# Create the time series plot
overall <- ggplot(data = kelp_mean, aes(x = year, y = kelp_mean)) +
  geom_line() +
  geom_ribbon(data = kelp_mean , aes(ymin = mean_kelp - sd_kelp, ymax = mean_kelp + sd_kelp), fill = "gray80", alpha = 0.5) +
  geom_point(data = subset_kelp, aes(color = color), size = 3) +
  geom_hline(data = site_stats , aes(yintercept = mean_kelp), linetype = "solid") + # add horizontal line
  geom_vline(data = site_stats , aes(xintercept = 2013), linetype = "dashed") + # add vertical line
  labs(x = "Year", y = "Kelp Mean", color = "") +
  scale_color_manual(values = c("Within 1 SD" = "black", "Above 1 SD" = "forestgreen", "Below 1 SD" = "purple")) +
  facet_wrap(~site, scales = "fixed", ncol=5) +
  theme_classic()+
  #add t-test results
  # Add p-values with color based on significance level
  geom_text(data = ttest_results, aes(label = paste("p-value:", signif(p_value, digits = 2)), 
                                      x = max(kelp_mean$year), y = max(kelp_mean$kelp_mean)), hjust = 1, vjust = 1)
  #geom_smooth(method = "lm", se = FALSE, color = "black", linetype= "dotted",formula = y ~ x)+
  #ggpmisc::stat_poly_eq(aes(label = paste("P-value:", signif(..p.value.., digits = 2))),
  #             formula = y ~ x, parse = TRUE)


overall



####split

# Calculate maximum values for each site

# Create the time series plot for resist_sites
panel_A <- ggplot(data = resist_sites, aes(x = year, y = kelp_mean)) +
  geom_line() +
  geom_ribbon(aes(ymin = mean_kelp - sd_kelp, ymax = mean_kelp + sd_kelp), fill = "gray80", alpha = 0.5) +
  geom_point(data = subset_kelp %>% filter(site %in% resist_sites$site), aes(color = color), size = 2) +
  geom_hline(data = site_stats %>% filter(site %in% resist_sites$site), aes(yintercept = mean_kelp), linetype = "solid") + # add horizontal line
  labs(x = "Year", y = "Kelp density \n(stipes per 60 m²)", color = "") +
  scale_color_manual(values = c("Within 1 SD" = "black", "Above 1 SD" = "forestgreen", "Below 1 SD" = "purple")) +
  facet_wrap(~reorder(site, -mean_kelp), scales = "fixed", nrow=1, ncol=5
             )+
  scale_y_continuous(breaks = c(0, 100,200,300), limits = c(0,350))+
  labs(tag = "A") +
  theme_classic() +
  theme(plot.margin = unit(c(0, 1, 0.5, 0), "lines"),
        aspect.ratio = 1.1) + my_theme +
  #add ttest p-val
  geom_text(data = resist_sites %>% filter(year == 2007), aes(label = paste("P[B,A]:", signif(p_value, digits = 2)), 
                                     x = 2013.5, y = 330, hjust = 1, vjust = 1), size=2)+
  # Heatwave
  annotate(geom="rect", xmin=2014, xmax=2016, ymin=-Inf, ymax=Inf, fill="indianred1", alpha=0.2) 

panel_A

# Create the time series plot for transition_sites
panel_B <- ggplot(data = transition_sites, aes(x = year, y = kelp_mean)) +
  geom_line() +
  geom_ribbon(aes(ymin = mean_kelp - sd_kelp, ymax = mean_kelp + sd_kelp), fill = "gray80", alpha = 0.5) +
  geom_point(data = subset_kelp %>% filter(site %in% transition_sites$site), aes(color = color), size = 2) +
  geom_hline(data = site_stats %>% filter(site %in% transition_sites$site), aes(yintercept = mean_kelp), linetype = "solid") + # add horizontal line
  labs(x = "Year", y = "Kelp density \n(stipes per 60 m²)", color = "") +
  scale_color_manual(values = c("Within 1 SD" = "black", "Above 1 SD" = "forestgreen", "Below 1 SD" = "purple")) +
  facet_wrap(~reorder(site, -mean_kelp), scales = "fixed", ncol=5) +
  scale_y_continuous(breaks = c(0, 100,200,300), limits = c(0,350))+
  labs(tag = "B") +
  theme_classic() +
  theme(plot.margin = unit(c(0, 1, 0.5, 0), "lines"),
        aspect.ratio = 1.1) + my_theme+
  #add ttest p-val
  geom_text(data = transition_sites %>% filter(year == 2007), aes(label = paste("P[B,A]:", signif(p_value, digits = 2)), 
                                                              x = 2014, y = 320, hjust = 1, vjust = 1), size=2)+
  # Heatwave
  annotate(geom="rect", xmin=2014, xmax=2016, ymin=-Inf, ymax=Inf, fill="indianred1", alpha=0.2) 
panel_B

# Combine the two plots into one figure with two panels
g <- ggpubr::ggarrange(panel_A, panel_B, ncol = 1, common.legend = TRUE, legend = "top",
                  align = "v", heights = c(1, 3.1)
                  ) # + theme(panel.spacing = unit(0.01, "lines"))


g

ggsave(g, filename=file.path(figdir, "FigX_starting_conditions_new3.png"), 
       width=7, height=8, units="in", dpi=600, bg="white")





