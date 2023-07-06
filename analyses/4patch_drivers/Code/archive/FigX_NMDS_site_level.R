#Communtiy change analyses
#Joshua G. Smith; jossmith@mbayaq.org

rm(list=ls())

librarian::shelf(tidyverse, here, vegan, ggplot2, cluster, ggforce)


################################################################################
#set directories and load data
basedir <- "/Volumes/seaotterdb$/kelp_recovery/"

stan_dat <- read.csv(file.path(basedir, "data/subtidal_monitoring/processed/kelp_stan_CC.csv")) 

figdir <- here::here("analyses","4patch_drivers","Figures")

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

#----------------process standardized data--------------------------------------

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
cent <- aggregate(cbind(NMDS1, NMDS2) ~ year + site, data = scrs, FUN = mean) 


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


################################################################################
#Step 2 - plot

# Theme
my_theme <-  theme(axis.text=element_text(size=6, color = "black"),
                   axis.text.y = element_text(angle = 90, hjust = 0.5, color = "black"),
                   axis.title=element_text(size=8, color = "black"),
                   plot.tag=element_blank(), #element_text(size=8),
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
                   strip.text = element_text(size = 6 ,face="bold", hjust=0, color = "black"),
)



# plot cent_newroids with ellipses around clusters
library(scales)


#cent_new$year_range <- scales::rescale(as.numeric(cent_new$year), to = c(0,1))
#year_breaks <- unique(cent_new$year)[seq(1, length(unique(cent_new$year)), 4)]

year_breaks <- as.numeric(c("2008","2011","2014","2017","2020"))

stan_trajectory <- ggplot(data = cent_new %>% 
                            mutate(site = str_to_title(gsub("_", " ", site)),
                                   site = str_replace(site, "Dc", "DC"),
                                   site = str_replace(site, "Uc", "UC")),
                          aes(NMDS1, NMDS2, color = year)) +
  geom_path(size = 1, lineend = "round") +
  scale_color_gradient(low = "green", high = "purple", name = "Year", 
                       breaks = year_breaks, labels = year_breaks,
                     guide = guide_colorbar(title.position = "top", title.hjust = 0.5, title.vjust = 1)) +
  stat_ellipse(aes(fill = factor(clust_new)), type = "norm", level = 0.8, geom = "polygon", alpha = 0.2) +
  scale_fill_manual(values = c("forestgreen", "purple")) +
  ggtitle("Kelp forest community structure") +
  facet_wrap(~site, scales = "fixed") +
  theme_bw() +
  labs(fill = "K-means \ncluster")+
  theme(plot.title = element_text(size = 16, face = "bold"),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        legend.position = "right") + my_theme + coord_equal()
stan_trajectory




ggsave(stan_trajectory, filename=file.path(figdir, "FigSX_site_level_trajectories.png"), bg = "white",
       width=7, height=7, units="in", dpi=600) 


























