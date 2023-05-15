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
  summarise(kelp_mean = mean(counts))



#join
swath_sub <- left_join(swath_sub1, kelp_mean, by=c("year", "site", "outbreak_period"))


################################################################################
#calculate % change of kelp and urchins at the transect level

df_mean_kelp <- swath_sub %>%
  group_by(site, species, zone, transect, outbreak_period) %>% 
  dplyr::summarise(mean_count = mean(counts),
            se_count = sd(counts, na.rm = TRUE)/sqrt(sum(!is.na(counts))))%>%
  pivot_wider(names_from = outbreak_period, values_from = c(mean_count, se_count)) %>%
  mutate_all(~ if_else(is.na(.), 0, .))

names(df_mean_kelp)

df_means <- df_mean_kelp %>%
              #add a constant
              mutate( mean_count_before = ifelse(mean_count_before == 0 , 0.0001,mean_count_before),
                perc_change = (mean_count_after-mean_count_before)/mean_count_before) 



#calculate mean % change of kelp and urchins at the site level
df_means2 <- df_means %>%
  group_by(site) %>%
  summarize(mean_macro = mean(perc_change[species == "macrocystis_pyrifera"], na.rm=TRUE),
            mean_strongy_perc = mean(perc_change[species == "strongylocentrotus_purpuratus"],na.rm=TRUE),
            mean_strongy = mean(mean_count_after[species == "strongylocentrotus_purpuratus"],na.rm=TRUE),
            se_macro = sd(perc_change[species == "macrocystis_pyrifera"]) / sqrt(sum(species == "macrocystis_pyrifera")),
            se_strongy_perc = sd(perc_change[species == "strongylocentrotus_purpuratus"]) / sqrt(sum(species == "strongylocentrotus_purpuratus")),
            se_strongy = sd(mean_count_after[species == "strongylocentrotus_purpuratus"]) / sqrt(sum(species == "strongylocentrotus_purpuratus"))
            )

################################################################################
#join with sim 

plot_dat <- left_join(df_means2, dist_long1, by="site")


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


#%change by density plot
p1 <- ggplot(plot_dat %>% mutate(site = str_to_title(gsub("_", " ", site)),
                           site = str_replace(site, "Dc", "DC"),
                           site = str_replace(site, "Uc", "UC")),
       aes(x = mean_strongy/60, y = mean_macro, label = site, color = mean_sim)) +
  geom_point(size=3) +
  geom_errorbar(aes(ymin = mean_macro - se_macro, ymax = mean_macro + se_macro), width = 0) +
  geom_errorbarh(aes(xmin = mean_strongy/60 - se_strongy/60, xmax = mean_strongy/60 + se_strongy/60), height = 0) +
  labs(color = "Community similarity \n (2007-2013 vs. 2014-2020)")+
  xlab("Mean purple sea urchin density 2014-2020 (per m²)") +
  ylab("Kelp stipe percent change \n(2014-2020 vs. 2007-2013)") +
  ggrepel::geom_label_repel(box.padding = 1, point.padding = 0, color = "black", size=2) +
  scale_color_gradient(low = "#B4AEE8", high = "#82B366")+
  theme_bw()+my_theme

ggsave(p1, filename=file.path(figdir, "Fig2_kelp_urchins_sim.png"), 
      width=7, height=6, units="in", dpi=600)







