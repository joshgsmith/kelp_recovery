
#Joshua G. Smith; jossmith@mbayaq.org

rm(list=ls())


######
#required packages
librarian::shelf(tidyverse, sf, zoo)

#set directories 
basedir <- "/Volumes/seaotterdb$/kelp_recovery/data/"
censusdir <- file.path(basedir,"census_data/annual_surveys/processed")
gisdir <- "/Volumes/seaotterdb$/kelp_recovery/data/gis_data"

#read census data
census_orig <- st_read(file.path(censusdir, "mpen_counts_1985_23_beta.geojson")) 



################################################################################


# Theme
base_theme <-  theme(axis.text=element_text(size=12, color = "black"),
                     axis.title=element_text(size=12,color = "black"),
                     plot.tag=element_text(size=7,color = "black"),
                     plot.title=element_text(size=10,color = "black", face = "bold"),
                     # Gridlines
                     panel.grid.major = element_blank(), 
                     panel.grid.minor = element_blank(),
                     panel.background = element_blank(), 
                     axis.line = element_line(colour = "black"),
                     # Legend
                     legend.key.size = unit(0.3, "cm"), 
                     legend.key = element_rect(fill=alpha('blue', 0)),
                     legend.spacing.y = unit(0.1, "cm"),  
                     legend.text=element_text(size=8,color = "black"),
                     #legend.title=element_(size=8,color = "black"),
                     #legend.key.height = unit(0.1, "cm"),
                     legend.background = element_rect(fill=alpha('blue', 0)),
                     #facets
                     strip.text = element_text(size=10, face = "bold",color = "black", hjust=0),
                     strip.background = element_blank())



census_join_summary <- census_orig %>%
  group_by(year) %>%
  summarize(total_n_indep = sum(n_indep), total_n_pup = sum(n_pup))

ggplot(census_join_summary, aes(x = year)) +
  geom_point(aes(y = total_n_indep, color = "Independents")) +  
  geom_point(aes(y = total_n_pup, color = "Pups")) +            
  geom_line(aes(y = total_n_indep, color = "Independents")) +  
  geom_line(aes(y = total_n_pup, color = "Pups")) +            
  labs(x = "Year", y = "Total Count") +
  scale_color_manual(values = c("Independents" = "blue", "Pups" = "red"),
                     guide = guide_legend()) +  
  ggtitle("Total Counts of Independents and Pups Over Time") +
  theme_minimal() +
  theme(legend.title = element_blank()) + base_theme


