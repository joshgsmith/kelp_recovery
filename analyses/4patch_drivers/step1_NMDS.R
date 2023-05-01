#Communtiy change analyses
#Joshua G. Smith; jossmith@mbayaq.org

rm(list=ls())

librarian::shelf(tidyverse, here, vegan)


################################################################################
#set directories and load data
basedir <- "/Volumes/seaotterdb$/kelp_recovery/"

stan_dat <- read.csv(file.path(basedir, "data/subtidal_monitoring/processed/kelp_stan_CC.csv")) 


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
#scores and plot

#Calculate centroids
scrs<- as.data.frame(scores(stan_ord, display="site"))

scrs <- cbind(as.data.frame(scrs), year=stan_group_vars$year, site = stan_group_vars$site) #to facilitate computing centroids, add group var
cent <- aggregate(cbind(NMDS1, NMDS2) ~ year + site, data = scrs, FUN = mean) #computes centroids by MHW and region




stan_trajectory<- ggplot(data=cent %>%
                           mutate(basin = ifelse(year < 2012, "before",ifelse(year > 2014, "after",NA)))
                    
                         , aes(NMDS1, NMDS2)) +
  #add year-to-year trajectory
  geom_segment(aes(x = NMDS1, y = NMDS2, xend = lead(NMDS1), yend = lead(NMDS2)),
               arrow = arrow(length = unit(0.25, "cm"), type = "closed"),
               size = 0.5, color = "black",
               data = cent, size=4) +
  #add confidence ellipses
  stat_ellipse(level=0.95, size=1.1, aes(color=basin))+
  ggtitle("Kelp forest community structure")+
  #scale_color_manual(values=c('#D81B60','#1E88E5'))+
  theme(strip.text = element_text(size=18, face="bold"))+
  theme(plot.title = element_text(size=16, face="bold"))+
  theme(text = element_text(20))+
  #labs(shape="Heatwave period",color='Site type')+
  #ggrepel::geom_label_repel(aes(label = year),
   #                         box.padding   = 1, 
    #                      point.padding = 0.1,
     #                      segment.color = 'grey',
      #                      max.overlaps=Inf,
       #                     size=4) +
  facet_wrap(~site, scales="free")+
  scale_color_manual(values=c("purple","forestgreen"))+
  theme_bw()

stan_trajectory




################################################################################
#cluster analysis to determine optimal grouping

# Create a list of data frames, one for each site
site_list <- split(stan_dat, stan_dat$site)

# Loop through each site and perform the cluster analysis and plot the dendrogram
for (site in names(site_list)) {
  # Subset the data to only include the current site
  site_data <- site_list[[site]]
  
  # Calculate the distance matrix using a suitable distance measure
  dist_matrix <- vegdist(site_data[, 10:ncol(site_data)], method = "bray")
  
  # Perform the cluster analysis
  cluster_result <- hclust(dist_matrix)
  
  # Cut the dendrogram to obtain the cluster assignments
  cluster_assignments <- cutree(cluster_result, k = length(unique(site_data$year)))
  
  # Add the cluster assignments as a new column to the original data
  site_data$clusters <- cluster_assignments
  
  # Create a vector of labels to be used in the plot
  label_vec <- site_data$year
  
  # Plot the dendrogram for this site, with custom labels
  plot(cluster_result, labels = label_vec, main = paste("Dendrogram for", site))
}
