#Communtiy change analyses
#Joshua G. Smith; jossmith@mbayaq.org

rm(list=ls())

librarian::shelf(tidyverse, here, vegan, ggplot2, cluster, ggforce, reshape2)


################################################################################
#set directories and load data
figdir <- here::here("analyses","4patch_drivers","Figures")
basedir <- here::here("analyses","4patch_drivers","Output")

#load multivariate data
load(file.path(basedir, "multivariate_data.Rdata"))

#load standardized data
stan_dat <- read.csv(file.path(basedir, "kelp_stan_CC.csv")) 

#load raw dat
fish_raw <- read.csv(file.path(basedir, "kelp_fish_counts_CC.csv")) 

upc_raw <- read.csv(file.path(basedir, "kelp_upc_cov_CC.csv")) 

swath_raw <- read.csv(file.path(basedir, "kelp_swath_counts_CC.csv")) %>%
  #remove kelps -- we will handle them as their own group
  dplyr::select(-macrocystis_pyrifera,
                -pterygophora_californica,
                -eisenia_arborea,
                -stephanocystis_osmundacea,
                -laminaria_setchellii,
                -nereocystis_luetkeana,
                -alaria_marginata,
                -costaria_costata,
                -laminaria_farlowii,
                -pleurophycus_gardneri)

kelp_raw <- read.csv(file.path(basedir, "kelp_swath_counts_CC.csv")) %>% 
  #extract kelps as their own group              
          dplyr::select(1:11, macrocystis_pyrifera,
                              pterygophora_californica,
                              eisenia_arborea,
                              stephanocystis_osmundacea,
                              laminaria_setchellii,
                              nereocystis_luetkeana,
                              alaria_marginata,
                              costaria_costata,
                              laminaria_farlowii,
                              pleurophycus_gardneri)


################################################################################
#process data

#replace any NAs with 0
stan_dat <- stan_dat %>% mutate(across(where(is.numeric), ~replace_na(., 0)))

fish_sum <- fish_raw %>% group_by(year, MHW, site) %>%
  dplyr::summarize(across(10:114, mean, na.rm = TRUE))

swath_sum <- swath_raw %>% group_by(year, MHW, site) %>%
  dplyr::summarize(across(9:57, mean, na.rm = TRUE))

upc_sum <- upc_raw %>% group_by(year, MHW, site) %>%
  dplyr::summarize(across(9:64, mean, na.rm = TRUE))

kelp_sum <- kelp_raw %>% group_by(year, MHW, site) %>%
  dplyr::summarize(across(9:18, mean, na.rm = TRUE))


################################################################################
#Step 2 - determine optimal centroid clustering

#indexing
species_names <- colnames(stan_ord_dat)
# Find the index of the species of interest
species_index <- match("strongylocentrotus_purpuratus", species_names)
species_scores <- scores(stan_ord, display = "sites")[,] 


scrs<- as.data.frame(scores(stan_ord, display="site"))
scrs <- cbind(as.data.frame(scrs), year=stan_group_vars$year) #to facilitate computing centroids, add group var
cent <- aggregate(cbind(NMDS1, NMDS2) ~ year, data = scrs, FUN = mean) #computes centroids by MHW and region

# determine optimal number of clusters using elbow method
wss <- (nrow(cent)-1)*sum(apply(cent[,2:3],2,var))
for(i in 1:13) wss[i] <- sum(kmeans(cent[,2:3],centers=i)$withinss)
plot(1:13, wss, type="b", xlab="Number of clusters", ylab="Within groups sum of squares")

# Create a data frame for plotting
df <- data.frame(NumClusters = 1:13, WSS = wss)

# Plot using ggplot2
elbow_plot <- ggplot(df, aes(x = NumClusters, y = WSS)) +
  geom_line(size=0.5) +
  geom_point(size=1) +
  labs(x = "Number of clusters", y = "Within groups sum of squares", size=2) +
  theme_bw() + theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank(),
                     axis.text = element_text(size = 6),
                     axis.title = element_text(size=8))


#ggsave(elbow_plot, filename=file.path(figdir, "FigSX_elbow_plot.png"), bg = "white",
 #      width=4, height=4, units="in", dpi=600) 


# choose number of clusters based on elbow in plot
k <- 2

# perform k-means clustering and add cluster label to data frame
set.seed(123)
cent$cluster <- kmeans(cent[,2:3], k)$cluster

################################################################################
#Step 3 - process diversity
#process swath data

swath_dat <- swath_sum %>% ungroup() %>% dplyr::select(4:ncol(.))
swath_groups <- swath_sum %>% ungroup() %>% dplyr::select(1:3)

swath_richness <- data.frame(S.obs = apply(swath_dat[,1:49]>0, 1, sum)) #count each unique species once
swath_evenness <- diversity(swath_dat)/log(specnumber(swath_dat))
swath_shannon <- diversity(swath_dat, index="shannon")
swath_simpson <- diversity(swath_dat, index="simpson")
swath_abund <- rowSums(swath_dat[,1:49])

swath_alphadiv <- cbind(swath_groups, swath_richness, swath_shannon, swath_simpson, swath_evenness, swath_abund)%>%
  mutate(MHW = str_to_sentence(MHW),
         MHW = factor(MHW, levels=c("Before","During","After")))


#process upc data
upc_dat <- upc_sum %>% ungroup() %>% dplyr::select(4:ncol(.))
upc_groups <- upc_sum %>% ungroup() %>% dplyr::select(1:3)

total_area_covered <- rowSums(upc_dat)

abundance <- upc_dat / total_area_covered * 166.67 #convert %cov to abundance, multiply by scalar, in this case 60m2 transect, 166.67 (which is 10,000 divided by 60)

upc_dat[,1:56] <- abundance

#convert perc cov to relative abundance

upc_richness <- data.frame(S.obs = apply(upc_dat[,1:56]>0, 1, sum))
upc_evenness <- diversity(upc_dat)/log(specnumber(upc_dat))
upc_shannon <- diversity(upc_dat, index="shannon")
upc_simpson <- diversity(upc_dat, index="simpson")
upc_abund <- rowSums(upc_dat[,1:56])

upc_alphadiv <- cbind(upc_groups, upc_richness, upc_shannon, upc_simpson, upc_evenness, upc_abund)%>%
  mutate(MHW = str_to_sentence(MHW),
         MHW = factor(MHW, levels=c("Before","During","After")))


#process fish data

fish_dat <- fish_sum %>% ungroup() %>% dplyr::select(4:ncol(.))
fish_groups <- fish_sum %>% ungroup() %>% dplyr::select(1:3)

fish_richness <- data.frame(S.obs = apply(fish_dat[,1:105]>0, 1, sum))
fish_evenness <- diversity(fish_dat)/log(specnumber(fish_dat))
fish_shannon <- diversity(fish_dat, index="shannon")
fish_simpson <- diversity(fish_dat, index="simpson")
fish_abund <- rowSums(fish_dat[,1:105])

fish_alphadiv <- cbind(fish_groups, fish_richness, fish_shannon, fish_simpson, fish_evenness, fish_abund)%>%
  mutate(MHW = str_to_sentence(MHW),
         MHW = factor(MHW, levels=c("Before","During","After")),
         fish_evenness = ifelse(fish_evenness == "NaN",NA,fish_evenness))


#process kelps

kelp_dat <- kelp_sum %>% ungroup() %>% dplyr::select(4:ncol(.))
kelp_groups <- kelp_sum %>% ungroup() %>% dplyr::select(1:3)

kelp_richness <- data.frame(S.obs = apply(kelp_dat[,1:10]>0, 1, sum)) #count each unique species once
kelp_evenness <- diversity(kelp_dat)/log(specnumber(kelp_dat))
kelp_shannon <- diversity(kelp_dat, index="shannon")
kelp_simpson <- diversity(kelp_dat, index="simpson")
kelp_abund <- rowSums(kelp_dat[,1:10])

kelp_alphadiv <- cbind(kelp_groups, kelp_richness, kelp_shannon, kelp_simpson, kelp_evenness, kelp_abund)%>%
  mutate(MHW = str_to_sentence(MHW),
         MHW = factor(MHW, levels=c("Before","During","After")))

################################################################################
#Step 4 determine spp that explain changes over time

fish_join <- cbind(fish_alphadiv, fish_dat) %>% rename(richness=S.obs) %>%
              #calculate annual mean
              group_by(year)%>%
              summarize(across(3:112,mean))

fish_rich <- fish_join %>%
  select(year, 7:111) %>%
  mutate(across(2:106, ~ifelse(. > 0, 1, 0))) 

# Create a data frame for species changes (added or lost)
species_changes <- data.frame(year = numeric(0), species_added = character(0), species_lost = character(0))

# Sort the data frame by year
fish_rich <- fish_rich %>%
  arrange(year)

# Iterate through years starting from the second year (2008)
for (i in 2:nrow(fish_rich)) {
  current_year <- fish_rich$year[i]
  previous_year <- fish_rich$year[i - 1]
  
  current_species <- fish_rich[i, 3:106]  # Exclude the year column from comparison
  previous_species <- fish_rich[i - 1, 3:106]
  
  species_added <- colnames(current_species)[current_species > previous_species]
  species_lost <- colnames(previous_species)[previous_species > current_species]
  
  species_changes <- rbind(species_changes, data.frame(year = current_year, species_added = paste(species_added, collapse = ", "), species_lost = paste(species_lost, collapse = ", ")))
}

# Print or work with species_changes to see which species were added or lost each year relative to the previous year
print(species_changes)


# Reshape the data to long format
fish_rich_long <- fish_rich %>%
  pivot_longer(cols = -year, names_to = "species", values_to = "presence")

# Filter out species that were never observed (presence == 0) in any year
fish_rich_filtered <- fish_rich_long %>%
  group_by(species) %>%
  filter(any(presence == 1))

# Calculate richness (number of species observed) for each year based on the filtered dataset
richness_data <- fish_rich_filtered %>%
  group_by(year) %>%
  summarize(richness = sum(presence))

# Filter out species that were observed in ALL years
fish_rich_reduced <- fish_rich_filtered %>%
  group_by(species) %>%
  filter(sum(presence) != n_distinct(year))

# Create a heatmap
g <- ggplot(data = fish_rich_reduced, aes(x = year, y = species, fill = factor(presence))) +
  geom_tile() +
  scale_fill_manual(values = c("0" = "white", "1" = "green"), labels = c("0" = "Absent", "1" = "Present")) +
  labs(title = "Species Presence/Absence Over Years", x = "Year", y = "Species") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))



# Create a colored line plot with a gradient color transition between years
ggplot(data = richness_data, aes(x = year, y = 0, color = richness)) +
  geom_line(size = 20, aes(group = 1)) +
  labs(title = "Richness Over Years", x = "Year", y = NULL) +
  scale_color_gradientn(
    colors = c("navyblue", "indianred"),
    values = scales::rescale(c(min(richness_data$year), max(richness_data$year))),
    guide = "none"
  ) +  # Adjust the gradient colors and values as needed
  theme_minimal()

################################################################################
#Step 5 - examine site cohesion and distance changes over time


stan_group_vars2 <- stan_group_vars %>% mutate(period = ifelse(year < 2013,"Before","After"),
                                               site_period = paste(site, period),
                                               year_site_period = paste(year, site, period)) %>% data.frame()

# Set row and column names in stan_max_distmat
stan_max_distmat2 <- usedist::dist_setNames(stan_max_distmat, stan_group_vars2$year_site_period)


#use betadisper to reduce vegdist to principal coords
disper_mat <- betadisper(stan_max_distmat2, type="centroid",
                         group = stan_group_vars2$year)

boxplot(disper_mat)

#create function to calculate distance between samples and centroid
betadistances <-
  function(x)
  {
    cnt <- x$centroids
    coord <- x$vectors
    pos <- which(x$eig >= 0)
    neg <- which(x$eig < 0)
    d <- apply(cnt[,pos], 1,
               function(z) rowSums(sweep(coord[,pos], 2, z)^2))
    if (length(neg))
      d <- d - apply(cnt[, neg], 1,
                     function(z) rowSums(sweep(coord[,neg], 2, z)^2))
    d <- as.data.frame(sqrt(d))
    cbind("group" = x$group, d)
  }

cen_distance <- betadistances(disper_mat) %>%
                  tibble::rownames_to_column(var = "sample_name") %>%
                  #reshape
                  mutate(centroid = word(sample_name, -1),
                         year = word(sample_name, 1),
                         site = word(sample_name, 2)) %>%
                  janitor::clean_names() %>%
                  dplyr::select(-group) %>%
                  dplyr::select(sample_name, centroid, year, site,everything()) %>%
                  #make longer
                  pivot_longer(cols = 5:52, names_to = "centroid_2", values_to = "distance") %>%
                   mutate(centroid_period = str_extract(centroid_2, "([^_]+)$"),
                          centroid_2 = str_remove(centroid_2, "_after"),
                          centroid_2 = str_remove(centroid_2, "_before"),
                          centroid_2 = toupper(centroid_2)) %>%
                  #filter distances to 'before' centroid
                  filter(centroid_period == "before") %>%
                  #make sure we compare within sites only by matching site and centroid
                  filter(site == centroid_2) %>%
                  #clean up
                  dplyr::select(year, site, distance) %>%
                  mutate(year = as.factor(year),
                         site = as.factor(site)) %>%
                  data.frame()

#check that it worked. 

nrow(cen_distance) == nrow(stan_group_vars2)

###the resulting output of cen_distance is the distance to the 'before' centroid
#for a given site over time. Centroids correspond to the site-level centroid
#of years 2007-2012. 

################################################################################
#Step 4 - plot


my_theme <-  theme(axis.text=element_text(size=10, color = "black"),
                   axis.text.y = element_text(angle = 90, hjust = 0.5, color = "black"),
                   axis.title=element_text(size=10, color = "black"),
                   plot.tag=element_text(size=8, color = "black"),
                   plot.title =element_text(size=8, face="bold", color = "black"),
                   # Gridlines 
                   panel.grid.major = element_blank(), 
                   panel.grid.minor = element_blank(),
                   panel.background = element_blank(), 
                   axis.line = element_line(colour = "black"),
                   # Legend
                   legend.key = element_blank(),
                   legend.background = element_rect(fill=alpha('blue', 0)),
                   legend.key.height = unit(1, "lines"), 
                   legend.text = element_text(size = 7, color = "black"),
                   legend.title = element_text(size = 8, color = "black"),
                   #legend.spacing.y = unit(0.75, "cm"),
                   #facets
                   strip.background = element_blank(),
                   strip.text = element_text(size = 6 ,face="bold", color ="black"),
)


#plot NMDS
stan_trajectory <- ggplot(data = cent %>% mutate(basin = ifelse(year < 2012, "before",ifelse(year > 2014, "after",NA))),
                          aes(NMDS1, NMDS2)) +
  geom_segment(aes(x = NMDS1, y = NMDS2, xend = lead(NMDS1), yend = lead(NMDS2)),
               arrow = arrow(length = unit(0.25, "cm"), type = "closed"),
               size = 0.5, color = "black",
               data = cent, size=4) +
  #geom_path(size = 1, lineend = "round") + #use geom_path to validate the segment
  stat_ellipse(aes(fill = factor(cluster)), type = "norm", level = 0.8, geom = "polygon", alpha = 0.2) +
  scale_fill_manual(values = c("forestgreen", "purple")) +
  #labs(shape="Heatwave period",color='Site type')+
  ggrepel::geom_label_repel(aes(label = year),
                            box.padding   = 1, 
                            point.padding = 0.1,
                            segment.color = 'grey',
                            max.overlaps=Inf,
                            size=2) +
  ggtitle("Kelp forest community structure") +
  labs(color = 'K-means \ncluster', fill = "K-means \ncluster", tag = "A")+
  theme_bw() + my_theme + theme(#plot.margin = margin(t=2,1,1,1, "lines"),
    legend.position = c(0.9, 0.2)) #+ theme(plot.margin = margin(40,0,40,0))

#stan_trajectory


ymax <- max(swath_alphadiv$swath_shannon, na.rm=TRUE)
shannon <- ggplot() +
  #kelp
  geom_point(aes(x = year, y = kelp_shannon, color = "Kelps", fill = "Kelps"), alpha = 0.1, size = 0.5, data = kelp_alphadiv %>% mutate(site = gsub("_", " ", site)), position = position_jitter(width = 0.1)) +
  geom_smooth(method = "auto", se = TRUE, size = 0.5, alpha = 0.4, aes(x = year, y = kelp_shannon, color = "Kelps", fill = "Kelps"), data = kelp_alphadiv %>% mutate(site = gsub("_", " ", site))) +
  #swath
  geom_point(aes(x = year, y = swath_shannon, color = "Mobile and \nconspicuous inverts", fill = "Mobile and \nconspicuous inverts"), alpha = 0.1, size = 0.5, data = swath_alphadiv %>% mutate(site = gsub("_", " ", site)), position = position_jitter(width = 0.1)) +
  geom_smooth(method = "auto", se = TRUE, size = 0.5, alpha = 0.4, aes(x = year, y = swath_shannon, color = "Mobile and \nconspicuous inverts", fill = "Mobile and \nconspicuous inverts"), data = swath_alphadiv %>% mutate(site = gsub("_", " ", site))) +
  #upc
  geom_point(aes(x = year + 0.1, y = upc_shannon, color = "Sessile inverts and \nmacroalgae", fill = "Sessile inverts and \nmacroalgae"), alpha = 0.1, size = 0.5, data = upc_alphadiv %>% mutate(site = gsub("_", " ", site)), position = position_jitter(width = 0.1)) +
  geom_smooth(method = "auto", se = TRUE, size = 0.5, alpha = 0.4, aes(x = year, y = upc_shannon, color = "Sessile inverts and \nmacroalgae", fill = "Sessile inverts and \nmacroalgae"), data = upc_alphadiv %>% mutate(site = gsub("_", " ", site))) +
  #fish
  geom_point(aes(x = year - 0.1, y = fish_shannon, color = "Fishes", fill = "Fishes"), alpha = 0.1, size = 0.5, data = fish_alphadiv %>% mutate(site = gsub("_", " ", site)), position = position_jitter(width = 0.1)) +
  geom_smooth(method = "auto", se = TRUE, size = 0.5, alpha = 0.4, aes(x = year, y = fish_shannon, color = "Fishes", fill = "Fishes"), data = fish_alphadiv %>% mutate(site = gsub("_", " ", site))) +
  #color
  scale_color_manual(values = c("Kelps" = "#009E73", "Mobile and \nconspicuous inverts" = "#E69F00", "Sessile inverts and \nmacroalgae" = "#7570B3", "Fishes" = "#CC79A7"), breaks = c("Sessile inverts and \nmacroalgae", "Mobile and \nconspicuous inverts", "Kelps", "Fishes")) +
  scale_fill_manual(values = c("Kelps" = "#009E73", "Mobile and \nconspicuous inverts" = "#E69F00", "Sessile inverts and \nmacroalgae" = "#7570B3", "Fishes" = "#CC79A7"), breaks = c("Sessile inverts and \nmacroalgae", "Mobile and \nconspicuous inverts", "Kelps", "Fishes")) +
  # Heatwave
  annotate(geom = "rect", xmin = 2014, xmax = 2016, ymin = -Inf, ymax = Inf, fill = "indianred1", alpha = 0.2) +
  theme_bw() + my_theme +
  labs(color = "Taxa", fill = "Taxa") +
  ylab("Shannon diversity") +
  xlab("Year") +
  ggtitle("Shannon diversity") +
  guides(color = guide_legend(keyheight = unit(1.1, "lines")), fill = guide_legend(keyheight = unit(1.1, "lines")))+
  labs(tag = "C") + theme(#axis.text.x = element_blank(), 
                          axis.title.x = element_blank(),
                          #plot.margin = margin(5, 0, 20, 3)
                          )


#shannon

ymax <- max(swath_alphadiv$swath_simpson, na.rm=TRUE)
simpson <- ggplot() +
  #kelp
  geom_point(aes(x = year, y = kelp_simpson, color = "Kelps", fill = "Kelps"), alpha = 0.1, size = 0.5, data = kelp_alphadiv %>% mutate(site = gsub("_", " ", site)), position = position_jitter(width = 0.1)) +
  geom_smooth(method = "auto", se = TRUE, size = 0.5, alpha = 0.4, aes(x = year, y = kelp_simpson, color = "Kelps", fill = "Kelps"), data = kelp_alphadiv %>% mutate(site = gsub("_", " ", site))) +
  #swath
  geom_point(aes(x = year, y = swath_simpson, color = "Mobile and \nconspicuous inverts", fill = "Mobile and \nconspicuous inverts"), alpha = 0.1, size = 0.5, data = swath_alphadiv %>% mutate(site = gsub("_", " ", site)), position = position_jitter(width = 0.1)) +
  geom_smooth(method = "auto", se = TRUE, size = 0.5, alpha = 0.4, aes(x = year, y = swath_simpson, color = "Mobile and \nconspicuous inverts", fill = "Mobile and \nconspicuous inverts"), data = swath_alphadiv %>% mutate(site = gsub("_", " ", site))) +
  #upc
  geom_point(aes(x = year + 0.1, y = upc_simpson, color = "Sessile inverts and \nmacroalgae", fill = "Sessile inverts and \nmacroalgae"), alpha = 0.1, size = 0.5, data = upc_alphadiv %>% mutate(site = gsub("_", " ", site)), position = position_jitter(width = 0.1)) +
  geom_smooth(method = "auto", se = TRUE, size = 0.5, alpha = 0.4, aes(x = year, y = upc_simpson, color = "Sessile inverts and \nmacroalgae", fill = "Sessile inverts and \nmacroalgae"), data = upc_alphadiv %>% mutate(site = gsub("_", " ", site))) +
  #fish
  geom_point(aes(x = year - 0.1, y = fish_simpson, color = "Fishes", fill = "Fishes"), alpha = 0.1, size = 0.5, data = fish_alphadiv %>% mutate(site = gsub("_", " ", site)), position = position_jitter(width = 0.1)) +
  geom_smooth(method = "auto", se = TRUE, size = 0.5, alpha = 0.4, aes(x = year, y = fish_simpson, color = "Fishes", fill = "Fishes"), data = fish_alphadiv %>% mutate(site = gsub("_", " ", site))) +
  #color
  scale_color_manual(values = c("Kelps" = "#009E73", "Mobile and \nconspicuous inverts" = "#E69F00", "Sessile inverts and \nmacroalgae" = "#7570B3", "Fishes" = "#CC79A7"), breaks = c("Sessile inverts and \nmacroalgae", "Mobile and \nconspicuous inverts", "Kelps", "Fishes")) +
  scale_fill_manual(values = c("Kelps" = "#009E73", "Mobile and \nconspicuous inverts" = "#E69F00", "Sessile inverts and \nmacroalgae" = "#7570B3", "Fishes" = "#CC79A7"), breaks = c("Sessile inverts and \nmacroalgae", "Mobile and \nconspicuous inverts", "Kelps", "Fishes")) +
  # Heatwave
  annotate(geom = "rect", xmin = 2014, xmax = 2016, ymin = -Inf, ymax = Inf, fill = "indianred1", alpha = 0.2) +
  theme_bw() + my_theme +
  labs(color = "Taxa", fill = "Taxa") +
  ylab("Simpson diversity") +
  xlab("Year") +
  ggtitle("Simpson diversity") +
  guides(color = guide_legend(keyheight = unit(1.1, "lines")), fill = guide_legend(keyheight = unit(1.1, "lines")))
#simpson


richness <- ggplot()+
  #kelp
  geom_point(aes(x = year, y=S.obs, color = "Kelps", fill="Kelps"), alpha = 0.1, size=0.5, data = kelp_alphadiv %>% mutate(site = gsub("_", " ", site)),
             position = position_jitter(width = 0.1)) +
  geom_smooth(method = "auto", se = TRUE, size = 0.5, alpha = 0.4, aes(x = year, y=S.obs, color = "Kelps", fill="Kelps"),data = kelp_alphadiv %>% mutate(site = gsub("_", " ", site))) +
  #swath
  geom_point(aes(x = year, y=S.obs, color = "Mobile and \nconspicuous inverts", fill="Mobile and \nconspicuous inverts"), alpha = 0.1, size=0.5, data = swath_alphadiv %>% mutate(site = gsub("_", " ", site)),
             position = position_jitter(width = 0.1)) +
  geom_smooth(method = "auto", se = TRUE, size = 0.5, alpha = 0.4, aes(x = year, y=S.obs, color = "Mobile and \nconspicuous inverts", fill="Mobile and \nconspicuous inverts"),data = swath_alphadiv %>% mutate(site = gsub("_", " ", site))) +
  #fish
  geom_point(aes(x = year-0.1, y=S.obs, color = "Fishes", fill="Fishes"), alpha = 0.1, size=0.5, data = fish_alphadiv %>% mutate(site = gsub("_", " ", site)), 
             position = position_jitter(width = 0.1)) +
  geom_smooth(method = "auto", se = TRUE, size = 0.5, alpha = 0.4, aes(x = year, y=S.obs, color = "Fishes", fill="Fishes"),data = fish_alphadiv %>% mutate(site = gsub("_", " ", site)) ) +
  #upc
  geom_point(aes(x = year+0.1, y=S.obs, color = "Sessile inverts and \nmacroalgae",fill="Sessile inverts and \nmacroalgae"), alpha = 0.1, size=0.5, data = upc_alphadiv %>% mutate(site = gsub("_", " ", site)), 
             position = position_jitter(width = 0.1)) +
  geom_smooth(method = "auto", se = TRUE, size = 0.5, alpha = 0.4, aes(x = year, y=S.obs, color = "Sessile inverts and \nmacroalgae", fill = "Sessile inverts and \nmacroalgae"),data = upc_alphadiv %>% mutate(site = gsub("_", " ", site))) +
  #color
  scale_color_manual(values = c("Kelps" = "#009E73", "Mobile and \nconspicuous inverts" = "#E69F00", "Sessile inverts and \nmacroalgae" = "#7570B3", "Fishes" = "#CC79A7"), breaks = c("Sessile inverts and \nmacroalgae", "Mobile and \nconspicuous inverts", "Kelps", "Fishes")) +
  scale_fill_manual(values = c("Kelps" = "#009E73", "Mobile and \nconspicuous inverts" = "#E69F00", "Sessile inverts and \nmacroalgae" = "#7570B3", "Fishes" = "#CC79A7"), breaks = c("Sessile inverts and \nmacroalgae", "Mobile and \nconspicuous inverts", "Kelps", "Fishes")) +
  # Heatwave
  annotate(geom="rect", xmin=2014, xmax=2016, ymin=-Inf, ymax=Inf, fill="indianred1", alpha=0.2) +
  theme_bw() + my_theme +
  labs(color = "Taxa", fill = "Taxa")+
  ylab("Taxonomic richness (n)")+
  xlab("Year")+
  ggtitle("Taxonomic richness")+
  guides(color = guide_legend(keyheight = unit(1.1, "lines")), fill = guide_legend(keyheight = unit(1.1, "lines")))+
  labs(tag = "D") + theme(#axis.text.x = element_blank(), 
                          axis.title.x = element_blank(),
                          #plot.margin = margin(5, 0, 20, 3)
                          ) 

#richness



evenness <- ggplot()+
  #kelp
  geom_point(aes(x = year, y=kelp_evenness, color = "Kelps", fill="Kelps"), alpha = 0.1, size=0.5, data = kelp_alphadiv %>% mutate(site = gsub("_", " ", site)),
             position = position_jitter(width = 0.1)) +
  geom_smooth(method = "auto", se = TRUE, size = 0.5, alpha = 0.4, aes(x = year, y=kelp_evenness, color = "Kelps", fill="Kelps"),data = kelp_alphadiv %>% mutate(site = gsub("_", " ", site))) +
  #swath
  geom_point(aes(x = year, y=swath_evenness, color = "Mobile and \nconspicuous inverts", fill="Mobile and \nconspicuous inverts"), alpha = 0.1, size=0.5, data = swath_alphadiv %>% mutate(site = gsub("_", " ", site)),
             position = position_jitter(width = 0.1)) +
  geom_smooth(method = "auto", se = TRUE, size = 0.5, alpha = 0.4, aes(x = year, y=swath_evenness, color = "Mobile and \nconspicuous inverts", fill="Mobile and \nconspicuous inverts"),data = swath_alphadiv %>% mutate(site = gsub("_", " ", site))) +
  #fish
  geom_point(aes(x = year-0.1, y=fish_evenness, color = "Fishes", fill="Fishes"), alpha = 0.1, size=0.5, data = fish_alphadiv %>% mutate(site = gsub("_", " ", site)), 
             position = position_jitter(width = 0.1)) +
  geom_smooth(method = "auto", se = TRUE, size = 0.5, alpha = 0.4, aes(x = year, y=fish_evenness, color = "Fishes", fill="Fishes"),data = fish_alphadiv %>% mutate(site = gsub("_", " ", site)) ) +
  #upc
  geom_point(aes(x = year+0.1, y=upc_evenness, color = "Sessile inverts and \nmacroalgae",fill="Sessile inverts and \nmacroalgae"), alpha = 0.1, size=0.5, data = upc_alphadiv %>% mutate(site = gsub("_", " ", site)), 
             position = position_jitter(width = 0.1)) +
  geom_smooth(method = "auto", se = TRUE, size = 0.5, alpha = 0.4, aes(x = year, y=upc_evenness, color = "Sessile inverts and \nmacroalgae", fill = "Sessile inverts and \nmacroalgae"),data = upc_alphadiv %>% mutate(site = gsub("_", " ", site))) +
  #color
  scale_color_manual(values = c("Kelps" = "#009E73", "Mobile and \nconspicuous inverts" = "#E69F00", "Sessile inverts and \nmacroalgae" = "#7570B3", "Fishes" = "#CC79A7"), breaks = c("Sessile inverts and \nmacroalgae", "Mobile and \nconspicuous inverts", "Kelps", "Fishes")) +
  scale_fill_manual(values = c("Kelps" = "#009E73", "Mobile and \nconspicuous inverts" = "#E69F00", "Sessile inverts and \nmacroalgae" = "#7570B3", "Fishes" = "#CC79A7"), breaks = c("Sessile inverts and \nmacroalgae", "Mobile and \nconspicuous inverts", "Kelps", "Fishes")) +
  # Heatwave
  annotate(geom="rect", xmin=2014, xmax=2016, ymin=-Inf, ymax=Inf, fill="indianred1", alpha=0.2) +
  theme_bw() + my_theme +
  labs(color = "Taxa", fill = "Taxa")+
  ylab("Taxonomic evenness (n)")+
  xlab("Year")+
  ggtitle("Taxonomic evenness")+
  scale_y_continuous(labels = scales::number_format(accuracy = 0.1)) +
  guides(color = guide_legend(keyheight = unit(1.1, "lines")), fill = guide_legend(keyheight = unit(0.5, "lines")))+
  labs(tag = "E") #+ theme(plot.margin = margin(-2, 0, 9, 3))

#evenness

# Combine ggplots with ggarrange()
combined_plot <- ggpubr::ggarrange(shannon, richness, evenness, ncol = 1, nrow=3, common.legend = TRUE, legend = "right") 

# Display the plot
#combined_plot <- combined_plot + theme(plot.margin = margin(5, 0, 5, 20, "pt")) + 
 # labs(tag = "B") + theme(plot.tag=element_text(size=8, color = "black", face = "plain")) 

#ggsave(combined_plot, filename=file.path(figdir, "Fig3_richness_diversity_evenness.png"), bg = "white",
 #      width=5.5, height=4.5, units="in", dpi=600) 


region_wide_plot <- ggpubr::ggarrange(stan_trajectory, combined_plot, ncol=2, widths = c(0.6,0.4)) 


#ggsave(region_wide_plot, filename=file.path(figdir, "Fig3_regional_metrics_new6.png"), bg = "white",
 #     width=5.5, height=8, units="in", dpi=600) 



cen_plot <- ggplot(cen_distance, aes(x = as.numeric(as.character(year)), y = distance, group = site)) +
  geom_line(aes(color = "Sites (n=24)"), alpha = 0.7) +
  geom_smooth(aes(color = "Median", group = 1), method = "loess", span = 0.4, alpha = 0.8) +
  labs(x = "Year", y = "Distance (Euclidean) to 2007-2012 centroid", tag = "B",
       title = "Temporal stability",
       color = "") +  # Set the legend title here
  scale_color_manual(values = c("Median" = "black","Sites (n=24)" = "lightblue")) +
  guides(color = guide_legend()) +  # Remove guide_legend(title = "Trajectory")
  annotate(geom = "rect", xmin = 2014, xmax = 2016, ymin = -Inf, ymax = Inf, fill = "indianred1", alpha = 0.2) +
  annotate(geom = "text", label = "MHW", x = 2017.5, y = 0.6, size = 3) +
  annotate("segment", x = 2016.5, y = 0.6, xend = 2015, yend = 0.6, arrow = arrow(type = "closed", length = unit(0.02, "npc"))) +
  geom_vline(xintercept = 2013, linetype = "dotted", size = 0.3) +
  annotate(geom = "text", label = "SSW", x = 2011, y = 0.6, size = 3) +
  annotate("segment", x = 2011.2, y = 0.585, xend = 2012.8, yend = 0.5, arrow = arrow(type = "closed", length = unit(0.02, "npc"))) +
  theme_bw() + my_theme +
  guides(color = guide_legend(override.aes = list(fill = NA)),
         linetype = guide_legend(override.aes = list(fill = NA))) +
  theme(legend.key = element_rect(fill = NA, color = "white"),
        legend.position = c(0.86, 0.16))  




##arrange and export

left <- ggpubr::ggarrange(stan_trajectory, cen_plot, ncol=1, heights = c(0.5,0.5)) 
left


final_plot <- ggpubr::ggarrange(left, combined_plot, ncol=2)

ggsave(final_plot, filename=file.path(figdir, "Fig3_regional_metrics_new8.png"), bg = "white",
   width=8, height=7.5, units="in", dpi=600) 











########work on plotting proportions





# Calculate the total density for each site at each year
site_year_total_density <- swath_long %>%
  mutate(transition_site = ifelse(site == "HOPKINS_UC" | site == "CANNERY_UC" |
                                    site == "SIREN" | site == "CANNERY_DC","no","yes"))%>%
  group_by(year, site, trophic_ecology, transition_site) %>%
  summarize(total_density = sum(density, na.rm=TRUE))

# Calculate the mean total density for each MHW level
mean_total_density <- site_year_total_density %>%
  group_by(year, trophic_ecology, transition_site) %>%
  summarize(mean_density = mean(total_density, na.rm=TRUE))

# Calculate the proportion using the mean total density
proportion_data <- mean_total_density %>%
  group_by(year, transition_site) %>%
  mutate(proportion = mean_density / sum(mean_density))

# Create the bar plot
ggplot(proportion_data, aes(x = year, y = proportion, fill = trophic_ecology)) +
  geom_bar(stat = "identity") +
  labs(title = "Total Proportion of Trophic Ecology Levels by MHW",
       x = "MHW",
       y = "Proportion") +
  facet_wrap(~transition_site)+
  scale_fill_brewer(palette = "Set3") + # You can choose a different color palette
  theme_minimal()














