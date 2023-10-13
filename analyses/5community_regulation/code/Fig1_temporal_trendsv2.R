
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

# Get rocky intertidal data
meta_dat <- readxl::read_xlsx(file.path(basedir,"intertidal_monitoring/raw/mytilus_cov_pis_den.xlsx"),sheet = 1)
mus_orig <- readxl::read_xlsx(file.path(basedir,"intertidal_monitoring/raw/mytilus_cov_pis_den.xlsx"),sheet = 2)
pis_orig <- readxl::read_xlsx(file.path(basedir,"intertidal_monitoring/raw/mytilus_cov_pis_den.xlsx"),sheet = 3)

#read foraging data
forage_orig <- read_csv(file.path(basedir,"/foraging_data/processed/foraging_data_2016_2023.csv"))

#read bathy
bathy_5m <- st_read(file.path(basedir, "gis_data/raw/bathymetry/contours_5m/contours_5m.shp")) %>% filter(CONTOUR == "-5")

#read state
ca_counties_orig <- st_read(file.path(basedir, "gis_data/raw/ca_county_boundaries/s7vc7n.shp")) 

# Get land
usa <- rnaturalearth::ne_states(country="United States of America", returnclass = "sf")
foreign <- rnaturalearth::ne_countries(country=c("Canada", "Mexico"), returnclass = "sf")


################################################################################
#Step 1 - filter rocky intertidal to focal study area

mus_build1 <- mus_orig %>% filter(latitude >= 36.47986 & latitude <= 36.64640) %>%
  mutate(ssw_period = ifelse(year <=2013, "before","after")) %>%
  #drop asilomar since there is only one year
  filter(!(marine_site_name %in% c("Asilomar","China Rocks"))) 

#check
unique(mus_build1$marine_site_name)


pis_build1 <- pis_orig %>% filter(latitude >= 36.47986 & latitude <= 36.64640) %>%
  mutate(ssw_period = ifelse(year <=2013, "before","after")) %>%
  #drop asilomar since there is only one year
  filter(!(marine_site_name %in% c("Asilomar","China Rocks"))) 

################################################################################
#Step 2 - merge

#otters
census_join_summary <- census_orig %>%
  group_by(year) %>%
  summarize(total_n_indep = sum(n_indep), total_n_pup = sum(n_pup)) %>% st_drop_geometry() %>%
  pivot_longer(cols = c(total_n_indep, total_n_pup), 
               names_to = "category", 
               values_to = "response")

#mussels
mus_build2 <- mus_build1 %>% dplyr::select(year, site = marine_site_name, response = percent_cover) %>%
  #calculate mean
  group_by(year) %>% summarize(response = mean(response)) %>% mutate(category = "mussels") 

#stars
pis_build2 <- pis_build1 %>% dplyr::select(year, site = marine_site_name, response = density_per_m2) %>%
  #calculate mean
  group_by(year) %>% summarize(response = mean(response)) %>%   mutate(category = "P. ochraceus") 

# Join the data
combined_data <- bind_rows(census_join_summary, mus_build2, pis_build2) %>%
                  mutate(category = factor(category)) %>%
                  filter(year >=2000)%>%
                  filter(!(category %in% c("total_n_pup"))) %>%
                  mutate(category = case_when(
                    category == "total_n_indep" ~ "Total independent otters",
                    category == "mussels" ~ "Mussels",
                    TRUE ~ category
                  ),
                  #set order
                  category = factor(category, levels = c("Total independent otters","Mussels","P. ochraceus"))
                  )
                  
################################################################################
#Step X - extract mussel dives

forage_build1 <- forage_orig %>% 
  #filter mussel dives
  filter(prey == "mus") %>%
  #select bout locations
  dplyr::select(year, month, day, bout, lat, long) %>% distinct() %>%
  filter(!is.na(lat))%>%
  st_as_sf(coords = c("long","lat"), crs = 4326)


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


# Build inset
g1_inset <-  ggplotGrob(
  ggplot() +
    # Plot land
    geom_sf(data=foreign, fill="grey80", color="white", lwd=0.3) +
    geom_sf(data=usa, fill="grey80", color="white", lwd=0.3) +
    # Plot box
    annotate("rect", xmin=-122.6, xmax=-121, ymin=36.2, ymax=37.1, color="black", fill=NA, lwd=0.6) +
    # Label regions
    #geom_text(data=region_labels, mapping=aes(y=lat_dd, label=region), x= -124.4, hjust=0, size=2) +
    # Labels
    labs(x="", y="") +
    # Crop
    coord_sf(xlim = c(-124.5, -117), ylim = c(32.5, 42)) +
    # Theme
    theme_bw() + base_theme +
    theme( plot.margin = unit(rep(0, 4), "null"),
           panel.margin = unit(rep(0, 4), "null"),
           panel.background = element_rect(fill='transparent'), #transparent panel bg
           # plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
           axis.ticks = element_blank(),
           axis.ticks.length = unit(0, "null"),
           axis.ticks.margin = unit(0, "null"),
           axis.text = element_blank(),
           axis.title=element_blank(),
           axis.text.y = element_blank())
)

# Create the "Monterey" text label
monterey_label <- data.frame(
  x = c(-121.9, -121.97), # x-coordinate for the upper right corner
  y = c(36.64, 36.54),  # y-coordinate for the upper right corner
  label = c("Monterey \nBay", "Carmel \nBay")
)

#rocky intertidal sites
rocky_sites <- pis_build1 %>% dplyr::select(monitoring_site = marine_site_name, latitude, longitude)%>%distinct()


g <- ggplot() +
  # Add landmarks
  geom_text(data = monterey_label, mapping = aes(x = x, y = y, label = label),
            size = 3, fontface = "bold") +
  # Add CA inset
  annotation_custom(grob = g1_inset, 
                    xmin = -122.01, 
                    xmax = -121.96,
                    ymin = 36.625) +
  #labs(fill = "Max kelp \nextent")+
  guides(
    fill = guide_legend(
      override.aes = list(size = 1),  # Adjust the legend point size as needed
      title = NULL  # Remove the legend title
    )
  )+
  #add foraging bouts for mussels
  ggnewscale::new_scale_fill()+
  geom_sf(
    data = forage_build1, #%>% filter(year == 2017),
    aes(fill = "Mussel forage \nbout"),
    size=0.1
  ) +
  labs(
    fill = NULL,  # Remove the title from the legend
    guide = guide_legend(title = NULL),  # Remove the legend title
    tag = "B"
  )+
  guides(
    fill = guide_legend(
      override.aes = list(size = 3),  # Adjust the legend point size as needed
      title = NULL  # Remove the legend title
    )
  )+
  #add land
  geom_sf(data = ca_counties_orig, fill = "gray", color = "gray80") +
  geom_sf(data = bathy_5m, fill = "black", color = "black") +
  #add rocky intertidal sites
  geom_point(
    data = rocky_sites,
    aes(x = longitude, y = ifelse(monitoring_site == "Point Lobos",latitude+.001,latitude)),
    shape = 24,  
    size = 3,
    fill = "orange"
  )+
  #add scale bar
  ggsn::scalebar(x.min = -121.99, x.max = -121.88, 
                 y.min = 36.519, y.max = 36.645,
                 #anchor=c(x=-124.7,y=41),
                 location="bottomright",
                 dist = 2, dist_unit = "km",
                 transform=TRUE, 
                 model = "WGS84",
                 st.dist=0.02,
                 st.size=2,
                 border.size=.5,
                 height=.02
  )+
  #add north arrow
  ggsn::north(x.min = -121.99, x.max = -121.88, 
              y.min = 36.519, y.max = 36.65,
              location = "topright", 
              scale = 0.05, 
              symbol = 10)+
  coord_sf(xlim = c(-121.99, -121.88), ylim = c(36.519, 36.645), crs = 4326)+
  labs(title = "", tag = "D")+
  theme_bw() + base_theme + theme(
    axis.title = element_blank(),
    legend.position=c(.8,.2),
    #legend.position = "top",  # Position the legend at the top
    #legend.justification = "right",  # Align the legend to the right
    #legend.box = "vertical",  # Box style for the legend
    #legend.margin = margin(t = -10, r = 10),  # Adjust the top and right margins for positioning
    axis.text = element_blank(),
    legend.background = element_rect(fill='transparent'),
    legend.key = element_rect(fill='transparent')
  )


g


p1 <- ggplot(combined_data %>% filter(category == "P. ochraceus"), aes(x = year, y = response)) +
  geom_point(alpha = 0.4, color = "#E377C2", show.legend = FALSE) +
  stat_smooth(geom = "line", size = 1, span = 0.6, color = "#E377C2", show.legend = FALSE) +
  stat_smooth(
    method = "loess",
    geom = "ribbon",
    alpha = 0.2,
    fill = "#E377C2",
    color = NA,
    span = 0.6,
    show.legend = FALSE
  ) +
  geom_vline(xintercept = 2013, linetype = "dotted", size = 0.6) +
  #SSW
  geom_vline(xintercept = 2013, linetype = "dotted", size=0.3)+
  annotate(geom="text", label="SSW", x=2011.5, y=0.25, size=2.5) +
  annotate("segment", x = 2011.8, y = 0.23, xend = 2012.7, yend = 0.2,
           arrow = arrow(type = "closed", length = unit(0.02, "npc")))+
  theme_minimal() +
  theme(legend.title = element_blank()) +
  theme_bw() +
  base_theme +
  scale_y_continuous(limits = c(0, NA), oob = scales::squish) +
  scale_x_continuous(limits = c(2000, 2023)) +
  scale_color_manual(values = "#E377C2") + 
  scale_fill_manual(values = "#E377C2") +
  labs(x = "Year", y = "Density (no. per m²)", title = "P. ochraceus", tag = "C") 

p1




p2 <- ggplot(combined_data %>% filter(category == "Mussels"), aes(x = year, y = response)) +
  geom_point(alpha = 0.4, color = "#FF7F0E", show.legend = FALSE) +
  stat_smooth(geom = "line", size = 1, span = 0.6, color = "#FF7F0E", show.legend = FALSE) +
  stat_smooth(
    method = "loess",
    geom = "ribbon",
    alpha = 0.2,
    fill = "#FF7F0E",
    color = NA,
    span = 0.6,
    show.legend = FALSE
  ) +
  geom_vline(xintercept = 2013, linetype = "dotted", size = 0.6) +
  labs(x = "Year", y = "Percent cover") +
  ggtitle("") +
  theme_minimal() +
  theme(legend.title = element_blank()) +
  theme_bw() +
  base_theme +
  scale_y_continuous(limits = c(0, NA), oob = scales::squish) +
  scale_x_continuous(limits = c(2000, 2023)) +
  scale_color_manual(values = "#FF7F0E") + 
  scale_fill_manual(values = "#FF7F0E")+
  labs(x = "", y = "Percent cover", title = "Mussels", tag = "B") +
  theme(axis.text.x = element_blank())

p2


p3 <- ggplot(combined_data %>% filter(category == "Total independent otters"), aes(x = year, y = response)) +
  geom_point(alpha = 0.4, color = "#2CA02C", show.legend = FALSE) +
  stat_smooth(geom = "line", size = 1, span = 0.6, color = "#2CA02C", show.legend = FALSE) +
  stat_smooth(
    method = "loess",
    geom = "ribbon",
    alpha = 0.2,
    fill = "#2CA02C",
    color = NA,
    span = 0.6,
    show.legend = FALSE
  ) +
  geom_vline(xintercept = 2013, linetype = "dotted", size = 0.6) +
  ggtitle("") +
  theme_minimal() +
  theme(legend.title = element_blank()) +
  theme_bw() +
  base_theme +
  scale_y_continuous(limits = c(300, NA), oob = scales::squish) +
  scale_x_continuous(limits = c(2000, 2023)) +
  scale_color_manual(values = "#2CA02C") + 
  scale_fill_manual(values = "#2CA02C")  +    
  labs(x = "", y = "Number of independents",title = "Sea otters", tag = "A") +
  theme(axis.text.x = element_blank())

p3


p <- gridExtra::grid.arrange(p3,p2,p1, ncol=1) 

p_final <- gridExtra::grid.arrange(p,g, ncol=2)




#save
ggsave(p_final, filename = file.path(figdir, "Fig1_timeseriesv2.png"), 
       width =7, height = 5.5, units = "in", dpi = 600)


