
#Joshua G. Smith; jossmith@mbayaq.org

rm(list=ls())


######
#required packages
librarian::shelf(tidyverse, sf, zoo)

#set directories 
basedir <- "/Volumes/seaotterdb$/kelp_recovery/data/"
censusdir <- file.path(basedir,"census_data/annual_surveys/processed")
gisdir <- "/Volumes/seaotterdb$/kelp_recovery/data/gis_data"
figdir <- here::here("analyses","5community_regulation","figures")

#read census data
census_orig <- st_read(file.path(censusdir, "mpen_counts_1985_23_beta.geojson")) 



################################################################################


# Theme
base_theme <-  theme(axis.text=element_text(size=8, color = "black"),
                     axis.title=element_text(size=8,color = "black"),
                     plot.tag=element_text(size=7,color = "black"),
                     plot.title=element_text(size=10,color = "black", face = "bold"),
                     # Gridlines
                     panel.grid.major = element_blank(), 
                     panel.grid.minor = element_blank(),
                     panel.background = element_blank(), 
                     axis.line = element_line(colour = "black"),
                     # Legend
                     legend.key.size = unit(0.3, "cm"), 
                     #legend.key = element_rect(fill = "white"), # Set it to transparent
                     legend.spacing.y = unit(0.1, "cm"),  
                     legend.text=element_text(size=8,color = "black"),
                     legend.title=element_blank(),
                     #legend.key.height = unit(0.1, "cm"),
                     #legend.background = element_rect(fill=alpha('blue', 0)),
                     #facets
                     strip.text = element_text(size=10, face = "bold",color = "black", hjust=0),
                     strip.background = element_blank())

# Rest of your code
census_join_summary <- census_orig %>%
  group_by(year) %>%
  summarize(total_n_indep = sum(n_indep), total_n_pup = sum(n_pup))

g <- ggplot(census_join_summary, aes(x = year)) +
  geom_point(aes(y = total_n_indep, color = "Independents")) +  
  geom_point(aes(y = total_n_pup, color = "Pups")) +            
  geom_line(aes(y = total_n_indep, color = "Independents")) +  
  geom_line(aes(y = total_n_pup, color = "Pups")) +   
  #SSW
  geom_vline(xintercept = 2013, linetype = "dotted", size=0.6)+
  annotate(geom="text", label="Sea star \nwasting", x=2010, y=640 , size=3) +
  annotate("segment", x = 2010, y = 600, xend = 2012.7, yend = 550,
           arrow = arrow(type = "closed", length = unit(0.02, "npc")))+
  labs(x = "Year", y = "Total Count") +
  scale_color_brewer(palette = "Dark2", guide = guide_legend()) +  # Use Dark2 palette
  ggtitle("") +
  theme_minimal() +
  theme(legend.title = element_blank()) + theme_bw() + base_theme


#save
ggsave(g, filename = file.path(figdir, "FigX_timeseries.png"), 
       width = 6, height = 3, units = "in", dpi = 600)





