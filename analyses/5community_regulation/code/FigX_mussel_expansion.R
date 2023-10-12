#Joshua G. Smith; jossmith@mbayaq.org

rm(list=ls())


######
#required packages
librarian::shelf(tidyverse)

#set directories 
basedir <- "/Volumes/seaotterdb$/kelp_recovery/data" 
gisdir <- file.path(basedir, "gis_data/processed")
figdir <- here::here("analyses","5community_regulation","figures")
output <- here::here("analyses","5community_regulation","output")

# Get rocky intertidal position data
meta_dat <- readxl::read_xlsx(file.path(basedir,"intertidal_monitoring/raw/mytilus_pisaster_position_data.xlsx"),sheet = 1)
mus_pos_raw <- readxl::read_xlsx(file.path(basedir,"intertidal_monitoring/raw/mytilus_pisaster_position_data.xlsx"),sheet = 2)
star_pos_raw <- readxl::read_xlsx(file.path(basedir,"intertidal_monitoring/raw/mytilus_pisaster_position_data.xlsx"),sheet = 3)


################################################################################
#Step 1 - filter to focal study area

mus_pos_build1 <- mus_pos_raw %>% filter(latitude >= 36.47986 & latitude <= 36.64640) %>%
                    mutate(ssw_period = ifelse(year <=2013, "before","after")) %>%
                    #drop asilomar since there is only one year
                    filter(!(intertidal_sitename %in% c("Asilomar","China Rocks")))

#check
unique(mus_pos_build1$intertidal_sitename)

################################################################################
#plot heatmap by site


tile_theme <-  theme(axis.text=element_text(size=6, color = "black"),
                     axis.title=element_text(size=7,color = "black"),
                     plot.tag=element_text(size=7,color = "black"),
                     plot.title=element_text(size=7,color = "black", face = "bold"),
                     # Gridlines
                     panel.grid.major = element_blank(), 
                     panel.grid.minor = element_blank(),
                     panel.background = element_blank(), 
                     axis.line = element_line(colour = "black"),
                     # Legend
                     legend.key.size = unit(0.3, "cm"), 
                     legend.key = element_rect(fill=alpha('blue', 0)),
                     legend.spacing.y = unit(0.1, "cm"),  
                     legend.text=element_text(size=5,color = "black"),
                     legend.title=element_text(size=6,color = "black"),
                     #legend.key.height = unit(0.1, "cm"),
                     legend.background = element_rect(fill=alpha('blue', 0)),
                     #facets
                     strip.text = element_text(size=6, face = "bold",color = "black", hjust=0),
                     strip.background = element_blank())

# Calculate the range of points for each site across all years and identify points unique to 'after'
#We use !any(ssw_period == "before") to check if there are no rows with 'ssw_period'
#equal to "before" within each group of 'intertidal_sitename,' 'transect,' and 'location.' 
#If there are no such rows, it means the points are unique to the 'after' period, and we label them as "Only After."

heatmap_data <- mus_pos_build1 %>%
  group_by(intertidal_sitename, transect, location) %>%
  summarize(
    location_range = max(location) - min(location),
    Period = factor(!any(ssw_period == "before"), levels = c(TRUE, FALSE), labels = c("Post-SSW", "Pre-SSW"))
  ) %>%
  arrange(intertidal_sitename, transect) %>%
  group_by(intertidal_sitename) %>%
  mutate(transect_seq = dense_rank(transect),
         intertidal_sitename = factor(intertidal_sitename,levels = c("Hopkins","Point Pinos","Stillwater","Point Lobos")),
         Period = factor(Period, levels = c("Pre-SSW","Post-SSW")))


# Create a heatmap with equal tile size and color points observed only in 'after' differently
p <- ggplot(heatmap_data, aes(x = transect_seq, y = location, fill = Period)) +
  geom_tile(width=1, height=1) +  # Set width and height to 1 for equal tile size
  scale_fill_manual(values = c("Pre-SSW" = "navyblue", "Post-SSW" = "indianred")) +
  labs(title = "",
       x = "Transect number",
       y = "Distance from high to \nlow intertidal (meters from baseline)",
       fill = NULL) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_y_reverse()+  # Flip the y-axis
  #ggforce::facet_col(vars(intertidal_sitename), scales = "free") +
  facet_grid(.~intertidal_sitename) +
  coord_fixed()+
  theme_bw()+tile_theme  + theme(legend.position = "top")+
  theme(plot.margin = margin(0, 0, 0, 0, "cm")) #reduce right margin
p


# Create arrow scheme
arrow_data <- data.frame(x = 0.5, y = 0.85)

scheme1 <- ggplot(arrow_data, aes(x = x, y = y)) +
  geom_segment(aes(xend = x, y = 0.85, yend = y - 0.15),
               arrow = arrow(type = "closed", length = unit(0.1, "inches")),
               size = 1) +
  labs(x = NULL, y = NULL) +  # Remove axis labels
  theme_void() +  # Remove axis lines and ticks
  annotate("text", x = 0.5, y = 0.92, label = "High \nintertidal", hjust = 0.5, vjust=5.5, size = 2) +
  annotate("text", x = 0.5, y = 0.68, label = "Low \nintertidal", hjust = 0.5, vjust=-0.2, size = 2)+
  theme(plot.margin = margin(0, 1, 0, 0, "cm")) #reduce left margin


# Merge
layout_matrix <- matrix(c(1,2), ncol=2)
g1_full <- gridExtra::grid.arrange(p, scheme1,
                                   layout_matrix=layout_matrix,
                                  widths=c(0.9, 0.1))
g1_full



ggsave(g1_full, filename = file.path(figdir, "FigX_mussel_expansion.png"), 
       width = 7, height = 4, units = "in", dpi = 600)


