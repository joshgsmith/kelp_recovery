#Communtiy change analyses
#Joshua G. Smith; jossmith@mbayaq.org

rm(list=ls())

librarian::shelf(tidyverse, here, vegan, ggplot2, cluster, ggforce)


################################################################################
#set directories and load data
figdir <- here::here("analyses","4patch_drivers","Figures")
basedir <- here::here("analyses","4patch_drivers","Output")

#load multivariate data
load(file.path(basedir, "multivariate_data.Rdata"))

#load standardized data
stan_dat <- read.csv(file.path(basedir, "kelp_stan_CC.csv")) 

#load raw count data
fish_raw <- read.csv(file.path(basedir, "kelp_fish_counts_CC.csv"))

upc_raw <- read.csv(file.path(basedir, "kelp_upc_cov_CC.csv"))

swath_raw <- read.csv(file.path(basedir, "kelp_swath_counts_CC.csv")) 


################################################################################
#process data

#replace any NAs with 0
stan_dat <- stan_dat %>% mutate(across(where(is.numeric), ~replace_na(., 0))) 

#aggregate to site-year level by taking mean counts across replicate transects
fish_sum <- fish_raw %>% group_by(year, MHW, site) %>%
  dplyr::summarize(across(10:114, mean, na.rm = TRUE))

swath_sum <- swath_raw %>% group_by(year, MHW, site) %>%
  dplyr::summarize(across(9:67, mean, na.rm = TRUE))

upc_sum <- upc_raw %>% group_by(year, MHW, site) %>%
  dplyr::summarize(across(9:64, mean, na.rm = TRUE))



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

swath_richness <- data.frame(S.obs = apply(swath_dat[,1:59]>0, 1, sum))
swath_evenness <- diversity(swath_dat)/log(specnumber(swath_dat))
swath_shannon <- diversity(swath_dat, index="shannon")
swath_simpson <- diversity(swath_dat, index="simpson")
swath_abund <- rowSums(swath_dat[,1:59])

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

################################################################################
#Step 4 - plot


my_theme <-  theme(axis.text=element_text(size=6, color = "black"),
                   axis.text.y = element_text(angle = 90, hjust = 0.5, color = "black"),
                   axis.title=element_text(size=8, color = "black"),
                   plot.tag=element_text(size=8, color = "black"),
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
  theme_bw() + my_theme

#stan_trajectory


ymax <- max(swath_alphadiv$swath_shannon, na.rm=TRUE)
shannon <- ggplot()+
  #swath
  geom_point(aes(x = year, y=swath_shannon, color = "Kelp and \nmobile inverts", fill="Kelp and \nmobile inverts"), alpha = 0.1, size=0.5, data = swath_alphadiv %>% mutate(site = gsub("_", " ", site)),
             position = position_jitter(width = 0.1)) +
  geom_smooth(method = "auto", se = TRUE, size = 0.5, alpha = 0.4, aes(x = year, y=swath_shannon, color = "Kelp and \nmobile inverts", fill="Kelp and \nmobile inverts"),data = swath_alphadiv %>% mutate(site = gsub("_", " ", site))) +
  #fish
  geom_point(aes(x = year-0.1, y=fish_shannon, color = "Fish", fill="Fish"), alpha = 0.1, size=0.5, data = fish_alphadiv %>% mutate(site = gsub("_", " ", site)), 
             position = position_jitter(width = 0.1)) +
  geom_smooth(method = "auto", se = TRUE, size = 0.5, alpha = 0.4, aes(x = year, y=fish_shannon, color = "Fish", fill="Fish"),data = fish_alphadiv %>% mutate(site = gsub("_", " ", site)) ) +
  #upc
  geom_point(aes(x = year+0.1, y=upc_shannon, color = "Macroalgae and \nsessile inverts",fill="Macroalgae and \nsessile inverts"), alpha = 0.1, size=0.5, data = upc_alphadiv %>% mutate(site = gsub("_", " ", site)), 
             position = position_jitter(width = 0.1)) +
  geom_smooth(method = "auto", se = TRUE, size = 0.5, alpha = 0.4, aes(x = year, y=upc_shannon, color = "Macroalgae and \nsessile inverts", fill = "Macroalgae and \nsessile inverts"),data = upc_alphadiv %>% mutate(site = gsub("_", " ", site))) +
  #color
  scale_color_manual(values = c("#009E73", "#E69F00","#CC79A7"), labels = c("Fish", "Kelp and \nmobile inverts","Macroalgae and \nsessile inverts"))+ 
  scale_fill_manual(values = c("#009E73", "#E69F00","#CC79A7"), labels = c("Fish", "Kelp and \nmobile inverts","Macroalgae and \nsessile inverts"))+ 
  # Heatwave
  annotate(geom="rect", xmin=2014, xmax=2016, ymin=-Inf, ymax=Inf, fill="indianred1", alpha=0.2) +
  theme_bw() + my_theme +
  labs(color = "Organism \ntype", fill = "Organism \ntype")+
  ylab("Shannon diversity")+
  xlab("Year")+
  ggtitle("Shannon diversity")+
  guides(color = guide_legend(keyheight = unit(1.1, "lines")), fill = guide_legend(keyheight = unit(1.1, "lines")))

#shannon

ymax <- max(swath_alphadiv$swath_simpson, na.rm=TRUE)
simpson <- ggplot()+
  #swath
  geom_point(aes(x = year, y=swath_simpson, color = "Kelp and \nmobile inverts", fill="Kelp and \nmobile inverts"), alpha = 0.1, size=0.5, data = swath_alphadiv %>% mutate(site = gsub("_", " ", site)),
             position = position_jitter(width = 0.1)) +
  geom_smooth(method = "auto", se = TRUE, size = 0.5, alpha = 0.4, aes(x = year, y=swath_simpson, color = "Kelp and \nmobile inverts", fill="Kelp and \nmobile inverts"),data = swath_alphadiv %>% mutate(site = gsub("_", " ", site))) +
  #fish
  geom_point(aes(x = year-0.1, y=fish_simpson, color = "Fish", fill="Fish"), alpha = 0.1, size=0.5, data = fish_alphadiv %>% mutate(site = gsub("_", " ", site)), 
             position = position_jitter(width = 0.1)) +
  geom_smooth(method = "auto", se = TRUE, size = 0.5, alpha = 0.4, aes(x = year, y=fish_simpson, color = "Fish", fill="Fish"),data = fish_alphadiv %>% mutate(site = gsub("_", " ", site)) ) +
  #upc
  geom_point(aes(x = year+0.1, y=upc_simpson, color = "Macroalgae and \nsessile inverts",fill="Macroalgae and \nsessile inverts"), alpha = 0.1, size=0.5, data = upc_alphadiv %>% mutate(site = gsub("_", " ", site)), 
             position = position_jitter(width = 0.1)) +
  geom_smooth(method = "auto", se = TRUE, size = 0.5, alpha = 0.4, aes(x = year, y=upc_simpson, color = "Macroalgae and \nsessile inverts", fill = "Macroalgae and \nsessile inverts"),data = upc_alphadiv %>% mutate(site = gsub("_", " ", site))) +
  #color
  scale_color_manual(values = c("#009E73", "#E69F00","#CC79A7"), labels = c("Fish", "Kelp and \nmobile inverts","Macroalgae and \nsessile inverts"))+ 
  scale_fill_manual(values = c("#009E73", "#E69F00","#CC79A7"), labels = c("Fish", "Kelp and \nmobile inverts","Macroalgae and \nsessile inverts"))+ 
  # Heatwave
  annotate(geom="rect", xmin=2014, xmax=2016, ymin=-Inf, ymax=Inf, fill="indianred1", alpha=0.2) +
  theme_bw() + my_theme +
  labs(color = "Organism \ntype", fill = "Organism \ntype")+
  ylab("Simpson diversity")+
  xlab("Year")+
  ggtitle("Simpson diversity")+
  guides(color = guide_legend(keyheight = unit(1.1, "lines")), fill = guide_legend(keyheight = unit(1.1, "lines")))

#simpson


richness <- ggplot()+
  #swath
  geom_point(aes(x = year, y=S.obs, color = "Kelp and \nmobile inverts", fill="Kelp and \nmobile inverts"), alpha = 0.1, size=0.5, data = swath_alphadiv %>% mutate(site = gsub("_", " ", site)),
             position = position_jitter(width = 0.1)) +
  geom_smooth(method = "auto", se = TRUE, size = 0.5, alpha = 0.4, aes(x = year, y=S.obs, color = "Kelp and \nmobile inverts", fill="Kelp and \nmobile inverts"),data = swath_alphadiv %>% mutate(site = gsub("_", " ", site))) +
  #fish
  geom_point(aes(x = year-0.1, y=S.obs, color = "Fish", fill="Fish"), alpha = 0.1, size=0.5, data = fish_alphadiv %>% mutate(site = gsub("_", " ", site)), 
             position = position_jitter(width = 0.1)) +
  geom_smooth(method = "auto", se = TRUE, size = 0.5, alpha = 0.4, aes(x = year, y=S.obs, color = "Fish", fill="Fish"),data = fish_alphadiv %>% mutate(site = gsub("_", " ", site)) ) +
  #upc
  geom_point(aes(x = year+0.1, y=S.obs, color = "Macroalgae and \nsessile inverts",fill="Macroalgae and \nsessile inverts"), alpha = 0.1, size=0.5, data = upc_alphadiv %>% mutate(site = gsub("_", " ", site)), 
             position = position_jitter(width = 0.1)) +
  geom_smooth(method = "auto", se = TRUE, size = 0.5, alpha = 0.4, aes(x = year, y=S.obs, color = "Macroalgae and \nsessile inverts", fill = "Macroalgae and \nsessile inverts"),data = upc_alphadiv %>% mutate(site = gsub("_", " ", site))) +
  #color
  scale_color_manual(values = c("#009E73", "#E69F00","#CC79A7"), labels = c("Fish", "Kelp and \nmobile inverts","Macroalgae and \nsessile inverts"))+ 
  scale_fill_manual(values = c("#009E73", "#E69F00","#CC79A7"), labels = c("Fish", "Kelp and \nmobile inverts","Macroalgae and \nsessile inverts"))+ 
  # Heatwave
  annotate(geom="rect", xmin=2014, xmax=2016, ymin=-Inf, ymax=Inf, fill="indianred1", alpha=0.2) +
  theme_bw() + my_theme +
  labs(color = "Organism \ntype", fill = "Organism \ntype")+
  ylab("Species richness (n)")+
  xlab("Year")+
  ggtitle("Species richness")+
  guides(color = guide_legend(keyheight = unit(1.1, "lines")), fill = guide_legend(keyheight = unit(1.1, "lines")))

#richness



evenness <- ggplot()+
  #swath
  geom_point(aes(x = year, y=swath_evenness, color = "Kelp and \nmobile inverts", fill="Kelp and \nmobile inverts"), alpha = 0.1, size=0.5, data = swath_alphadiv %>% mutate(site = gsub("_", " ", site)),
             position = position_jitter(width = 0.1)) +
  geom_smooth(method = "auto", se = TRUE, size = 0.5, alpha = 0.4, aes(x = year, y=swath_evenness, color = "Kelp and \nmobile inverts", fill="Kelp and \nmobile inverts"),data = swath_alphadiv %>% mutate(site = gsub("_", " ", site))) +
  #fish
  geom_point(aes(x = year-0.1, y=fish_evenness, color = "Fish", fill="Fish"), alpha = 0.1, size=0.5, data = fish_alphadiv %>% mutate(site = gsub("_", " ", site)), 
             position = position_jitter(width = 0.1)) +
  geom_smooth(method = "auto", se = TRUE, size = 0.5, alpha = 0.4, aes(x = year, y=fish_evenness, color = "Fish", fill="Fish"),data = fish_alphadiv %>% mutate(site = gsub("_", " ", site)) ) +
  #upc
  geom_point(aes(x = year+0.1, y=upc_evenness, color = "Macroalgae and \nsessile inverts",fill="Macroalgae and \nsessile inverts"), alpha = 0.1, size=0.5, data = upc_alphadiv %>% mutate(site = gsub("_", " ", site)), 
             position = position_jitter(width = 0.1)) +
  geom_smooth(method = "auto", se = TRUE, size = 0.5, alpha = 0.4, aes(x = year, y=upc_evenness, color = "Macroalgae and \nsessile inverts", fill = "Macroalgae and \nsessile inverts"),data = upc_alphadiv %>% mutate(site = gsub("_", " ", site))) +
  #color
  scale_color_manual(values = c("#009E73", "#E69F00","#CC79A7"), labels = c("Fish", "Kelp and \nmobile inverts","Macroalgae and \nsessile inverts"))+ 
  scale_fill_manual(values = c("#009E73", "#E69F00","#CC79A7"), labels = c("Fish", "Kelp and \nmobile inverts","Macroalgae and \nsessile inverts"))+ 
  # Heatwave
  annotate(geom="rect", xmin=2014, xmax=2016, ymin=-Inf, ymax=Inf, fill="indianred1", alpha=0.2) +
  theme_bw() + my_theme +
  labs(color = "Organism \ntype", fill = "Organism \ntype")+
  ylab("Species evenness")+
  xlab("Year")+
  ggtitle("Species evenness")+
  guides(color = guide_legend(keyheight = unit(1.1, "lines")), fill = guide_legend(keyheight = unit(1.1, "lines")))
#evenness


# Combine ggplots with ggarrange()
combined_plot <- ggpubr::ggarrange(shannon, simpson, richness, evenness, ncol = 2, nrow=2, common.legend = TRUE, legend = "right") 

# Display the plot
combined_plot <- combined_plot + theme(plot.margin = margin(5, 0, 5, 20, "pt")) + 
  labs(tag = "B") + theme(plot.tag=element_text(size=8, color = "black", face = "plain")) 



region_wide_plot <- ggpubr::ggarrange(stan_trajectory, combined_plot, ncol=1)
region_wide_plot 

ggsave(region_wide_plot, filename=file.path(figdir, "Fig3_regional_metrics_new4.png"), bg = "white",
      width=5.5, height=8, units="in", dpi=600) 





