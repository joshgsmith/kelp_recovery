#Communtiy change analyses
#Joshua G. Smith; jossmith@mbayaq.org

rm(list=ls())

librarian::shelf(tidyverse, here, vegan, ggplot2, cluster, gganimate)


################################################################################
#set directories and load data
basedir <- "/Volumes/seaotterdb$/kelp_recovery/"
figdir <- "/Users/jossmith/kelp_recovery/analyses/4patch_drivers/Figures/gganimate"

stan_dat <- read.csv(file.path(basedir, "data/subtidal_monitoring/processed/kelp_stan_CC.csv")) 

#data available from https://opc.dataone.org/view/doi:10.25494/P6/MLPA_kelpforest.4

#processing script available at https://github.com/joshgsmith/kelp_recovery/tree/main/data/subtidal_monitoring

################################################################################
#pre-analysis processing

#replace any NAs with 0
stan_dat <- stan_dat %>% mutate(across(where(is.numeric), ~replace_na(., 0))) %>%
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
#prepare data for ordination

#define group vars
stan_group_vars <- stan_dat %>% dplyr::select(1:9)

#define data for ordination
stan_ord_dat <- stan_dat %>% dplyr::select(10:ncol(.))

#standardize to max 
stan_rel <- decostand(stan_ord_dat, method = "hellinger")

#generate a BC mat with stan dat
stan_max_distmat <- vegdist(stan_rel, method = "bray", na.rm = T)

#generate a BC mat with ord dat
stan_untransformed_distmat <- vegdist(stan_ord_dat, method = "bray", na.rm = T)


################################################################################
#ordinate data

set.seed(1985)
num_cores = 8

#ordinate stan dat
stan_ord <- metaMDS(stan_max_distmat, distance = "bray", parallel = num_cores, trymax=300)
stan_untrans_ord <- metaMDS(stan_untransformed_distmat, distance = "bray", parallel = num_cores, trymax=300)



################################################################################
#Step 2 - determine optimal centroid clustering

scrs<- as.data.frame(scores(stan_ord, display="site"))
scrs <- cbind(as.data.frame(scrs), year=stan_group_vars$year, site = stan_group_vars$site) #to facilitate computing centroids, add group var
cent <- aggregate(cbind(NMDS1, NMDS2) ~ year + site, data = scrs, FUN = mean) #computes centroids by year and site


#Sort cent by site and year to ensure clusters are assigned correctly
cent <- cent[order(cent$site, cent$year),] %>% mutate(site = as.factor(site))

# choose number of clusters based on elbow in plot
k <- 2

# perform k-means clustering and add cluster label to data frame
set.seed(123)
clusters <- lapply(split(cent, cent$site), function(subset) {
  kmeans(subset[,3:4], k)
})

cent$cluster <- unlist(lapply(clusters, function(x) x$cluster))

#check clusters. Smaller years should always be in cluster 1. 
wrong_clusters <- cent %>%
  group_by(site) %>%
  mutate(small_year = min(year)) %>%
  filter(small_year == year & cluster == 2) 

cent_new <- cent %>%
  mutate(clust_new = ifelse(site == "BLUEFISH_UC" |
                              site == "BUTTERFLY_DC" |
                              site == "CANNERY_DC"|
                              site == "HOPKINS_DC" |
                              site == "LOVERS_DC"|
                              site == "LOVERS_UC" |
                              site == "MACABEE_DC" |
                              site == "MONASTERY_DC"|
                              site == "MONASTERY_UC"|
                              site =="PESCADERO_UC"|
                              site == "SIREN" |
                              site == "STILLWATER_DC"|
                              site == "STILLWATER_UC"|
                              site == "WESTON_UC"
                            ,ifelse(cluster == 2,1,2),cluster))


cent_mean <- cent %>% group_by(year) %>% summarize(NMDS1 = mean(NMDS1),
                                                   NMDS2 = mean(NMDS2))

################################################################################
#plot


# Theme
my_theme <-  theme(axis.text=element_text(size=2),
                   axis.text.y = element_text(angle = 90, hjust = 0.5),
                   axis.title=element_text(size=2),
                   axis.ticks = element_line(size = 0.1),
                   plot.tag=element_text(size=2),
                   plot.title =element_text(size=2, margin = margin(0,0,7,0)),
                   plot.subtitle =element_text(size=3, face = "bold", margin = margin(0,0,-8,0), hjust = 0.1),
                   # Gridlines 
                   panel.grid.major = element_blank(), 
                   panel.grid.minor = element_blank(),
                   panel.background = element_blank(), 
                   #axis.line = element_line(colour = "black"),
                   # Legend
                   legend.key = element_blank(),
                   legend.background = element_rect(fill=alpha('blue', 0)),
                   legend.key.height = unit(1, "lines"), 
                   legend.text = element_text(size = 1),
                   legend.title = element_text(size = 2),
                   #legend.spacing.y = unit(0.75, "cm"),
                   #facets
                   strip.background = element_blank(),
                   strip.text = element_text(size = 1 ,face="bold"),
)


plot <- ggplot(cent_new, aes(x = NMDS1, y = NMDS2)) +
  geom_point(color = ifelse(cent_new$year %in% c(2013, 2014, 2015, 2016), "red", "#666666"), size=0.3) +
  geom_point(data=cent_mean, aes(x = NMDS1, y = NMDS2),color = "purple", size=3, alpha=0.6) +
  #transition_reveal(year)+
  transition_time(year) +
  ease_aes('linear')+
  theme_minimal()+
  labs(title = "Kelp forest community structure \nMonterey and Carmel Bay, CA")+
  labs(subtitle = "Year: {frame_time}")+
  #annotate("text", x = -0.4, y = 0.5, label = paste("Year:", "{frame_time}"))+
  shadow_mark(alpha = 0.2, size = 0.5)+
  annotate("text", x = 0.1, y = -0.395, label = paste("Created by: Joshua G. Smith"), 
           size=0.6, hjust = 0) +
  annotate("text", x = 0.1, y = -0.47, label = paste("Description: Relative abundance of algae, invertebrates, \nand fishes at 24 PISCO sites. Each point depicts a single \nlong-term site moving through multivariate space with the \ncentroid included as the purple bubble. Marine heatwave \nyears (2014-2016) are colored red."), 
           size=0.6, hjust = 0, vjust = 0.6) +
  annotate("text", x = 0.1, y = -0.58, label = paste("Data source: 10.25494/P6/MLPA_KELPFOREST.4"), 
           size=0.5, hjust = 0, fontface = "italic") +
  theme_bw()+my_theme 
plot

################################################################################
#export as .gif

animation <- animate(plot, renderer = gifski_renderer(), fps=7,
                     width = 1080,
                     height = 1080,
                     res = 600)


# Save the animation as a GIF file
anim_save(animation, filename=file.path(figdir, "comm_change_anim_credit.gif"))


################################################################################
#export as .mp4

# Set up the animation
animation <- animate(plot, renderer = av_renderer(), fps = 7,
                     width = 1080, height = 1080, res = 600)

# Save the animation as a video file
anim_save(animation,filename=file.path(figdir, "comm_change_anim.mp4"))


