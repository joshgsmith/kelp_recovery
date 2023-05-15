#Community data processing
#Joshua G. Smith
#May 1, 2023

rm(list=ls())
librarian::shelf(tidyverse, here, vegan, reshape2)


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

#community structure data for dist mat
stan_dat <- read.csv(file.path(basedir, "data/subtidal_monitoring/processed/kelp_stan_CC.csv")) %>%
  mutate(across(where(is.numeric), ~replace_na(., 0))) %>%
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
#generate distance matrix

set.seed(1985)

#define group vars
stan_group_vars <- stan_dat %>% dplyr::select(year, site) %>%
  #define urchin outbreak period
  mutate(period = ifelse(year > 2013, "after","before"),
         identifier = paste(year, site, period)) %>%
  dplyr::select(identifier)

#define data for matrix
stan_ord_dat <- stan_dat %>% dplyr::select(10:ncol(.))

#create matrix
stan_distmat <- as.matrix(vegdist(stan_ord_dat, method = "bray", na.rm = T))

#join group vars
dist_dat <- cbind(stan_group_vars, stan_distmat)

#create header names to match square matrix
colnames(dist_dat)[2:ncol(dist_dat)] <- dist_dat[,1]

#convert to three column format
dist_long <- setNames(melt(dist_dat), c('site_year_period1', 'site_year_period2', 'dissim')) 

dist_long1 <- dist_long %>% 
  #bust apart the grouping var
  mutate(year_1 = word(site_year_period1,1),
         year_2 = word(site_year_period2,1),
         period_1 = word(site_year_period1,-1),
         period_2 = word(site_year_period2,-1),
         site_1 = word(site_year_period1,-2),
         site_2 = word(site_year_period2,-2)) %>%
  #filter site-by-site comparisons only
  filter(site_1==site_2,
         period_1 == "after",
         period_2 == "before"
  ) %>%
  mutate(sim = 1-dissim)%>%
  #drop duplicates
  distinct()%>%
  #caclulate mean sim per site
  group_by(site_1)%>%
  summarize(mean_sim = mean(sim)) %>%
  rename(site = site_1)



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
  summarise(kelp_mean = mean(counts, na.rm=TRUE))



#join
swath_sub <- left_join(swath_sub1, kelp_mean, by=c("year", "site", "outbreak_period"))


################################################################################
#calculate % change of kelp and urchins at the transect level

df_mean_kelp <- swath_sub %>%
  group_by(site, species, outbreak_period) %>% 
  dplyr::summarise(mean_count = mean(counts, na.rm=TRUE),
                   se_count = sd(counts, na.rm = TRUE)/sqrt(sum(!is.na(counts))))%>%
  pivot_wider(names_from = outbreak_period, values_from = c(mean_count, se_count)) %>%
  mutate_all(~ if_else(is.na(.), 0, .))


################################################################################
#join with sim 

plot_dat <- left_join(df_mean_kelp, dist_long1, by="site")


################################################################################
#plot

# Theme
my_theme <-  theme(axis.text=element_text(size=6),
                   axis.text.y = element_text(angle = 90, hjust = 0.5),
                   axis.title=element_text(size=8),
                   plot.tag=element_blank(), #element_text(size=8),
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



# Subset the data by species
macro <- subset(plot_dat, species == "macrocystis_pyrifera") %>%
              #fix names
        mutate(site = str_to_title(gsub("_", " ", site)),
         site = str_replace(site, "Dc", "DC"),
         site = str_replace(site, "Uc", "UC"))
strong <- subset(plot_dat, species == "strongylocentrotus_purpuratus")%>%
  #fix names
  mutate(site = str_to_title(gsub("_", " ", site)),
         site = str_replace(site, "Dc", "DC"),
         site = str_replace(site, "Uc", "UC"))


# Create the dumbbell plot with arrows
A <- ggplot(data = macro, aes(x = mean_count_after, y = mean_count_before)) +
  geom_segment(aes(x = strong$mean_count_after, xend = strong$mean_count_after, y = macro$mean_count_before, yend = macro$mean_count_after, color = mean_sim), size = 1, arrow = arrow(length = unit(0.25, "cm"))) +
  geom_point(aes(x = strong$mean_count_after, y = macro$mean_count_before, color = mean_sim), size = 3) +
  xlab("Mean purple sea urchin density 2014-2020 (per 60 m²)") +
  ylab("Kelp stipe density (per 60 m²)") +
  labs(color = "Community similarity \n (2007-2013 vs. 2014-2020)")+
  #scale_x_log10("Mean Count After (Strongylocentrotus Purpuratus)") +
  scale_color_gradient(name = "Mean Sim") +
  scale_color_viridis_c() +
  theme_classic()+
  ggrepel::geom_label_repel(data = macro, aes(x = strong$mean_count_after, y = macro$mean_count_before, label = site),size=2, box.padding = 1, force = 50,min.segment.length = 1)+
  my_theme



# Create the dumbbell plot with arrows and log scale for strongylocentrotus_purpuratus
B <- ggplot(data = macro, aes(x = mean_count_after/60, y = mean_count_before/60)) +
  geom_segment(aes(x = strong$mean_count_before/60, xend = strong$mean_count_after/60, y = macro$mean_count_before/60, yend = macro$mean_count_after/60, color = mean_sim), size = 1, arrow = arrow(length = unit(0.25, "cm"))) +
  geom_point(aes(x = strong$mean_count_before/60, y = macro$mean_count_before/60, color = mean_sim), size = 3) +
  xlab("Purple sea urchin density (per m²)") +
  ylab("Kelp stipe density (per m²)") +
  labs(color = "Community similarity \n (2007-2013 vs. 2014-2020)")+
  #scale_x_log10("Purple sea urchin density (log no. per 60 m²)") +
  scale_color_viridis_c() +
  theme_classic()+
  ggrepel::geom_label_repel(data = macro, aes(x = strong$mean_count_before/60, y = macro$mean_count_before/60, label = site), size=2, box.padding = 1, force = 20,min.segment.length = 1)+
  my_theme

B

C <- ggpubr::ggarrange(A, B, ncol=1)
C


ggsave(B, filename=file.path(figdir, "Fig3_dumbbell.png"), 
       width=7, height=5, bg="white", units="in", dpi=600)





