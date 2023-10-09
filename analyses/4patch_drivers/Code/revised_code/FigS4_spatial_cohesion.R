#Communtiy change analyses
#Joshua G. Smith; jossmith@mbayaq.org

rm(list=ls())

require(librarian)

librarian::shelf(tidyverse, here, vegan, ggplot2, cluster, ggforce, data.table,
                 patchwork, scales)


################################################################################
#set directories and load data
basedir <- "/Volumes/seaotterdb$/kelp_recovery/"
outdir <- here::here("analyses","4patch_drivers","Output")
figdir <- here::here("analyses","4patch_drivers","Figures")

load(file.path(outdir,"multivariate_data.Rdata"))


################################################################################
#determine optimal centroid clustering

scrs<- as.data.frame(scores(stan_ord, display="site"))
scrs <- cbind(as.data.frame(scrs), year=stan_group_vars$year, site = stan_group_vars$site) #to facilitate computing centroids, add group var
cent <- aggregate(cbind(NMDS1, NMDS2) ~ year + site, data = scrs, FUN = mean) 


#Sort cent by site and year to ensure clusters are assigned correctly
cent <- cent[order(cent$site, cent$year),] %>% mutate(site = as.factor(site))


################################################################################
#test for cohesion

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

# Fill the lower triangle of the symmetric matrix
sym_matrix[lower.tri(sym_matrix)] <- t(sym_matrix)[lower.tri(sym_matrix)]

# Convert the symmetric matrix to a data frame
df <- reshape2::melt(sym_matrix)

df <- df %>% mutate(
  site_order_1 = case_when(
    Var1 == "CANNERY_UC" ~ 1,
    Var1 == "CANNERY_DC" ~ 2,
    Var1 == "MACABEE_UC" ~ 3,
    Var1 == "MACABEE_DC" ~ 4,
    Var1 == "HOPKINS_UC" ~ 5,
    Var1 == "HOPKINS_DC" ~ 6,
    Var1 == "LOVERS_UC" ~ 7,
    Var1 == "LOVERS_DC" ~ 8,
    Var1 == "SIREN" ~ 9,
    Var1 == "OTTER_PT_UC" ~ 10,
    Var1 == "OTTER_PT_DC" ~ 11,
    Var1 == "LONE_TREE" ~ 12,
    Var1 == "PESCADERO_UC" ~ 13,
    Var1 == "PESCADERO_DC" ~ 14,
    Var1 == "STILLWATER_UC" ~ 15,
    Var1 == "STILLWATER_DC" ~ 16,
    Var1 == "BUTTERFLY_UC" ~ 17,
    Var1 == "BUTTERFLY_DC" ~ 18,
    Var1 == "MONASTERY_UC" ~ 19,
    Var1 == "MONASTERY_DC" ~ 20,
    Var1 == "BLUEFISH_UC" ~ 21,
    Var1 == "BLUEFISH_DC" ~ 22,
    Var1 == "WESTON_UC" ~ 23,
    Var1 == "WESTON_DC" ~ 24,
    TRUE ~ NA
  ),
  site_order_2 = case_when(
    Var2 == "CANNERY_UC" ~ 1,
    Var2 == "CANNERY_DC" ~ 2,
    Var2 == "MACABEE_UC" ~ 3,
    Var2 == "MACABEE_DC" ~ 4,
    Var2 == "HOPKINS_UC" ~ 5,
    Var2 == "HOPKINS_DC" ~ 6,
    Var2 == "LOVERS_UC" ~ 7,
    Var2 == "LOVERS_DC" ~ 8,
    Var2 == "SIREN" ~ 9,
    Var2 == "OTTER_PT_UC" ~ 10,
    Var2 == "OTTER_PT_DC" ~ 11,
    Var2 == "LONE_TREE" ~ 12,
    Var2 == "PESCADERO_UC" ~ 13,
    Var2 == "PESCADERO_DC" ~ 14,
    Var2 == "STILLWATER_UC" ~ 15,
    Var2 == "STILLWATER_DC" ~ 16,
    Var2 == "BUTTERFLY_UC" ~ 17,
    Var2 == "BUTTERFLY_DC" ~ 18,
    Var2 == "MONASTERY_UC" ~ 19,
    Var2 == "MONASTERY_DC" ~ 20,
    Var2 == "BLUEFISH_UC" ~ 21,
    Var2 == "BLUEFISH_DC" ~ 22,
    Var2 == "WESTON_UC" ~ 23,
    Var2 == "WESTON_DC" ~ 24,
    TRUE ~ NA
))

# Create the heatmap using ggplot2
heatmap <- ggplot(df %>% mutate(Var1 = str_to_title(gsub("_", " ", Var1)),
                                Var1 = str_replace(Var1, "Dc", "DC"),
                                Var1 = str_replace(Var1, "Uc", "UC"),
                                Var2 = str_to_title(gsub("_", " ", Var2)),
                                Var2 = str_replace(Var2, "Dc", "DC"),
                                Var2 = str_replace(Var2, "Uc", "UC"),
                                Var1 = case_when(
                                  Var1 == "Cannery UC" ~ "Cannery UC*",
                                  Var1 == "Cannery DC" ~ "Cannery DC*",
                                  Var1 == "Hopkins UC" ~ "Hopkins UC*",
                                  Var1 == "Siren" ~ "Siren*",
                                  TRUE ~ Var1),
                                Var2 = case_when(
                                  Var2 == "Cannery UC" ~ "Cannery UC*",
                                  Var2 == "Cannery DC" ~ "Cannery DC*",
                                  Var2 == "Hopkins UC" ~ "Hopkins UC*",
                                  Var2 == "Siren" ~ "Siren*",
                                  TRUE ~ Var2)), aes(x = fct_reorder(Var1, 
                                                site_order_1),
                                                y = fct_reorder(Var2, site_order_2),
                                                fill = value))+
 # geom_tile() +
  geom_tile(data = . %>% filter(Var1 != Var2), na.rm = TRUE) +  # Filter out the lower triangle values
  scale_fill_gradient(low = "navyblue", high = "indianred", na.value = "white") +
  labs(x = "Site", y = "Site", fill = "Cohesion", tag = "B") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, # face = ifelse(df$Var1 %in% c("HOPKINS_UC", "CANNERY_UC", "SIREN", "CANNERY_DC"), "bold", "plain")
                               ),
    axis.text.y = element_text(hjust = 0.5, #face = ifelse(df$Var1 %in% c("HOPKINS_UC", "CANNERY_UC", "SIREN", "CANNERY_DC"), "bold", "plain")
                               ),
    legend.position = "right",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()) + theme_1 + theme(aspect.ratio = 0.55,
                                                          plot.margin = margin(3,5,5,5))


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
                            ,ifelse(cluster == 2,1,2),cluster),
         site_order = case_when(
           site == "CANNERY_UC" ~ 1,
           site == "CANNERY_DC" ~ 2,
           site == "MACABEE_UC" ~ 3,
           site == "MACABEE_DC" ~ 4,
           site == "HOPKINS_UC" ~ 5,
           site == "HOPKINS_DC" ~ 6,
           site == "LOVERS_UC" ~ 7,
           site == "LOVERS_DC" ~ 8,
           site == "SIREN" ~ 9,
           site == "OTTER_PT_UC" ~ 10,
           site == "OTTER_PT_DC" ~ 11,
           site == "LONE_TREE" ~ 12,
           site == "PESCADERO_UC" ~ 13,
           site == "PESCADERO_DC" ~ 14,
           site == "STILLWATER_UC" ~ 15,
           site == "STILLWATER_DC" ~ 16,
           site == "BUTTERFLY_UC" ~ 17,
           site == "BUTTERFLY_DC" ~ 18,
           site == "MONASTERY_UC" ~ 19,
           site == "MONASTERY_DC" ~ 20,
           site == "BLUEFISH_UC" ~ 21,
           site == "BLUEFISH_DC" ~ 22,
           site == "WESTON_UC" ~ 23,
           site == "WESTON_DC" ~ 24,
           TRUE ~ NA))

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
                   strip.text = element_text(size = 6 ,face="plain", hjust=0, color = "black"),
)


year_breaks <- as.numeric(c("2008","2011","2014","2017","2020"))

library(scales)

stan_trajectory <- ggplot(data = cent_new %>% 
                            mutate(site = str_to_title(gsub("_", " ", site)),
                                   site = str_replace(site, "Dc", "DC"),
                                   site = str_replace(site, "Uc", "UC"),
                                   site = case_when(
                                     site == "Cannery UC" ~ "Cannery UC*",
                                     site == "Cannery DC" ~ "Cannery DC*",
                                     site == "Hopkins UC" ~ "Hopkins UC*",
                                     site == "Siren" ~ "Siren*",
                                     TRUE ~ site)),
                          aes(NMDS1, NMDS2, color = year))+
  geom_path(size = 1, lineend = "round") +
  scale_color_gradient(low = "green", high = "purple", name = "Year", 
                       breaks = year_breaks, labels = year_breaks,
                       guide = guide_colorbar(title.position = "top", title.hjust = 0.5, title.vjust = 1)) +
  stat_ellipse(aes(fill = factor(clust_new)), type = "norm", level = 0.8, geom = "polygon", alpha = 0.2) +
  scale_fill_manual(values = c("forestgreen", "purple")) +
  facet_wrap(~reorder(site, site_order)) +
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
  coord_equal() + theme(plot.margin = margin(5,5,1,5),
                        plot.tag = element_text(vjust = -35)) 
stan_trajectory



# Arrange the plots side by side with equal heights and specified widths
arranged_plots <- stan_trajectory / plot_spacer() / heatmap +
  plot_layout(widths = c(4,-1.1,2),
              heights = c(2,-0.39,1)
              #heights = c(2,-0.1,1)
              )

arranged_plots


# Save the plot
ggsave(filename = file.path(figdir, "FigX_cohesion_new4.png"), plot = arranged_plots, 
       width = 6, height = 9.3, bg = "white", units = "in", dpi = 600)

#change height back to 9 potentially


