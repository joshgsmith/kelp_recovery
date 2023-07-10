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


################################################################################
#ordinate data

set.seed(1985)
num_cores = 8

#ordinate stan dat
stan_ord <- metaMDS(stan_max_distmat, distance = "bray", parallel = num_cores, trymax=300)


################################################################################
#compute centroids for years within each site. 

#Step 2 - determine optimal centroid clustering

scrs<- as.data.frame(scores(stan_ord, display="site"))
scrs <- cbind(as.data.frame(scrs), year=stan_group_vars$year, site = stan_group_vars$site) #to facilitate computing centroids, add group var
cent <- aggregate(cbind(NMDS1, NMDS2) ~ year + site, data = scrs, FUN = mean) 


#Sort cent by site and year to ensure clusters are assigned correctly
cent <- cent[order(cent$site, cent$year),] %>% mutate(site = as.factor(site))


################################################################################
#test for cohesion
library(vegan)  # for Procrustes analysis

# Empty lists to store results
cor_site <- list()
pvalue_site <- list()

# Create a vector of all unique years from 2007 to 2020
all_years <- 2007:2020

# Perform pairwise Procrustes analysis for each site pair
unique_sites <- unique(cent$site)
num_sites <- length(unique_sites)

# Interpolate missing years for each site
interpolated_cent <- list()
for (i in 1:num_sites) {
  site_i <- cent[cent$site == unique_sites[i], c("year", "NMDS1", "NMDS2")]
  
  # Identify missing years
  missing_years <- setdiff(all_years, site_i$year)
  
  if (length(missing_years) > 0) {
    # Interpolate missing years
    interpolated_pos <- data.frame(year = missing_years, NMDS1 = NA, NMDS2 = NA)
    site_i <- rbind(site_i, interpolated_pos)
    site_i <- site_i[order(site_i$year), ]  # Sort by year
    
    # Linearly interpolate NMDS1 and NMDS2 positions
    site_i$NMDS1 <- approx(site_i$year, site_i$NMDS1, xout = all_years)$y
    site_i$NMDS2 <- approx(site_i$year, site_i$NMDS2, xout = all_years)$y
  }
  
  # Preserve site name in interpolated_cent
  interpolated_cent[[unique_sites[i]]] <- site_i
}

# Create a data.table to store the results
library(data.table)
result_table <- data.table(site1 = character(),
                           site2 = character(),
                           correlation = numeric(),
                           p_value = numeric())

# Perform pairwise Procrustes analysis
for (i in 1:(num_sites - 1)) {
  site_i <- interpolated_cent[[unique_sites[i]]]
  
  for (j in (i+1):num_sites) {
    site_j <- interpolated_cent[[unique_sites[j]]]
    
    # Perform Procrustes analysis between site_i and site_j
    proc_pair <- protest(site_i[, c("NMDS1", "NMDS2")], site_j[, c("NMDS1", "NMDS2")])
    
    # Extract the correlation and p-value from the protest output
    correlation <- proc_pair$t0
    p_value <- proc_pair$signif
    
    # Add the results to the data.table
    result_table <- rbind(result_table, data.table(site1 = unique_sites[i],
                                                   site2 = unique_sites[j],
                                                   correlation = correlation,
                                                   p_value = p_value))
  }
}

################################################################################
#Plot

# Theme
theme_1 <-  theme(axis.text=element_text(size=6, color = "black"),
                   #axis.text.y = element_text(angle = 90, hjust = 0.5, color = "black"),
                    axis.text.x = element_text(angle = 45, hjust = 1, color = "black"),
                   axis.title=element_text(size=8, color = "black"),
                   plot.tag=element_text(size = 8, color = "black"), #element_text(size=8),
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

# Create a symmetric matrix with site names as row and column names
sym_matrix <- with(result_table, tapply(correlation, list(site1, site2), mean, na.rm = TRUE))
sites <- unique(c(result_table$site1, result_table$site2))
sym_matrix <- sym_matrix[sites, sites]

# Convert the symmetric matrix to a data frame
df <- reshape2::melt(sym_matrix)

# Create the heatmap using ggplot2
heatmap <- ggplot(df %>% mutate(Var1 = str_to_title(gsub("_", " ", Var1)),
                                Var1 = str_replace(Var1, "Dc", "DC"),
                                Var1 = str_replace(Var1, "Uc", "UC"),
                                Var2 = str_to_title(gsub("_", " ", Var2)),
                                Var2 = str_replace(Var2, "Dc", "DC"),
                                Var2 = str_replace(Var2, "Uc", "UC")), aes(x = Var1, y = Var2, fill = value)) +
  geom_tile() +
  scale_fill_gradient(low = "navyblue", high = "indianred", na.value = "white") +
  labs(x = "Site", y = "Site", fill = "Cohesion", tag = "B") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1),
    axis.text.y = element_text(hjust = 0.5),
    legend.position = "right",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()) + theme_1 + theme(aspect.ratio = 0.55,
                                                          plot.margin = margin(0,5,5,5))
  

# Display the heatmap
print(heatmap)


####plots trajectories

#scrs<- as.data.frame(scores(stan_ord, display="site"))
#scrs <- cbind(as.data.frame(scrs), year=stan_group_vars$year, site = stan_group_vars$site) #to facilitate computing centroids, add group var
#cent <- aggregate(cbind(NMDS1, NMDS2) ~ year + site, data = scrs, FUN = mean) 


#Sort cent by site and year to ensure clusters are assigned correctly
#cent <- cent[order(cent$site, cent$year),] %>% mutate(site = as.factor(site))

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


# Theme
my_theme <-  theme(axis.text=element_text(size=6, color = "black"),
                   axis.text.y = element_text(angle = 90, hjust = 0.5, color = "black"),
                   axis.title=element_text(size=8, color = "black"),
                   plot.tag=element_text(size= 8, color = "black"), #element_text(size=8),
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



year_breaks <- as.numeric(c("2008","2011","2014","2017","2020"))

library(scales)

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
  facet_wrap(~site) +
  theme_bw() +
  labs(fill = "K-means \ncluster", tag = "A") +
  theme(plot.title = element_text(size = 16, face = "bold"),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        legend.position = "right") + 
  my_theme + 
  theme(axis.title.y = element_text(vjust = 0)) +
  scale_x_continuous(labels = function(x) round(x, 1), breaks = pretty_breaks(n = 4)) +  # Round x-axis ticks to one decimal place and limit to 4 ticks
  scale_y_continuous(labels = function(y) round(y, 1), breaks = pretty_breaks(n = 4))+    # Round y-axis ticks to one decimal place and limit to 4 ticks
  coord_equal() + theme(plot.margin = margin(5,5,0,5),
                        plot.tag = element_text(vjust = -35)) 
stan_trajectory




library(patchwork)

# Arrange the plots side by side with equal heights and specified widths
arranged_plots <- stan_trajectory / plot_spacer() / heatmap +
  plot_layout(widths = c(4,-1.1,2),
              heights = c(2,-0.39,1)
              #heights = c(2,-0.1,1)
              )

arranged_plots


# Save the plot
ggsave(filename = file.path(figdir, "FigX_cohesion_new3.png"), plot = arranged_plots, 
       width = 6, height = 9, bg = "white", units = "in", dpi = 600)

#change height back to 9 potentially


