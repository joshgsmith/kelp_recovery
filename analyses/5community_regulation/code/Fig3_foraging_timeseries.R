

#Joshua G. Smith; jossmith@mbayaq.org

rm(list=ls())


######
#required packages

librarian::shelf(tidyverse, sf, raster, terra, janitor)

basedir <- "/Volumes/seaotterdb$/kelp_recovery/data"
figdir <- here::here("analyses","5community_regulation","figures")

#read landsat dat
landsat_orig <- st_read(file.path(basedir,"kelp_landsat/processed/monterey_peninsula/landsat_mpen_1984_2023_points_withNAs.shp"))

#read foraging data
forage_orig <- read_csv(file.path(basedir,"/foraging_data/processed/foraging_data_2016_2023.csv"))

#read sofa output
mass_class <- readxl::read_excel(file.path(basedir,"/sofa_data/raw/MCResults_Periods_10-11-23_12h.xlsx"), sheet = 4, skip=1) %>% clean_names()
dietcomp_class <- readxl::read_excel(file.path(basedir,"/sofa_data/raw/MCResults_Periods_10-11-23_12h.xlsx"), sheet = 5, skip = 1) %>% clean_names()
kcal_class <- readxl::read_excel(file.path(basedir,"/sofa_data/raw/MCResults_Periods_10-11-23_12h.xlsx"), sheet = 6, skip=1) %>% clean_names()

#read bathy
bathy_5m <- st_read(file.path(basedir, "gis_data/raw/bathymetry/contours_5m/contours_5m.shp")) %>% filter(CONTOUR == "-5")

#read state
ca_counties_orig <- st_read(file.path(basedir, "gis_data/raw/ca_county_boundaries/s7vc7n.shp")) 

# Get land
usa <- rnaturalearth::ne_states(country="United States of America", returnclass = "sf")
foreign <- rnaturalearth::ne_countries(country=c("Canada", "Mexico"), returnclass = "sf")


################################################################################
#determine max kelp extent as proxy for shallow subtidal


# transform landsat data to Teale Albers
rast_build1 <- st_transform(landsat_orig, crs = 3310) %>% 
  mutate(biomass = ifelse(biomass == 0, NA, biomass))  %>% filter (year == 2022) 


#select MPEN as focal region 
plot_dat <- rast_build1 %>% filter (latitude >= 36.510140 &
                                      latitude <= 36.670574) %>% st_transform(crs=3310)

#define 0 cover as historical kelp footprint
na_dat <- landsat_orig %>% filter (latitude >= 36.510140 &
                                     latitude <= 36.670574, 
                                   biomass == 0)

#transform landsat data to Teale Albers
t_dat <- st_transform(plot_dat, crs=3310) %>% 
  mutate(biomass = ifelse(biomass==0,NA,biomass))

kelp_historic <- st_transform(na_dat, crs=3310) # %>% mutate(biomass = ifelse(biomass==0,1,NA))

#create grid
r <- rast(t_dat, res=30)  # Builds a blank raster of given dimensions and resolution  
vr <- rasterize(t_dat, r,"biomass", resolution = 30) 

kelp_na <- rasterize(kelp_historic, r,"biomass", resolution = 30) 


################################################################################
#Step 2 - extract mussel dives

forage_build1 <- forage_orig %>% 
                    #filter mussel dives
                    filter(prey == "mus") %>%
                    #select bout locations
                    dplyr::select(year, month, day, bout, lat, long) %>% distinct() %>%
                    filter(!is.na(lat))%>%
                    st_as_sf(coords = c("long","lat"), crs = 4326)

################################################################################
#Step 3 - process sofa

ncol(mass_class)
ncol(dietcomp_class)
ncol(kcal_class)

# Add a column to each dataset indicating the source
mass_class1 <- mass_class %>% mutate(source = "mass_gain")
dietcomp_class1 <- dietcomp_class %>% mutate(source = "dietcomp")
kcal_class1 <- kcal_class %>% mutate(source = "kcal")

# Combine the datasets into one
sofa_build1 <- bind_rows(mass_class1, dietcomp_class1, kcal_class1)

################################################################################
#Step 4 - plot


base_theme <-  theme(axis.text=element_text(size=7, color = "black"),
                     axis.text.x = element_blank(),
                     axis.text.y = element_blank(),
                     axis.title=element_text(size=8,color = "black"),
                     legend.text=element_text(size=7,color = "black"),
                     legend.title=element_text(size=8,color = "black"),
                     plot.tag=element_text(size=10,color = "black"),
                     # Gridlines
                     panel.grid.major = element_blank(), 
                     panel.grid.minor = element_blank(),
                     panel.background = element_blank(), 
                     axis.line = element_line(colour = "black"),
                     # Legend
                     legend.key = element_rect(fill=alpha('blue', 0)),
                     legend.background = element_rect(fill=alpha('blue', 0)),
                     #facets
                     strip.text = element_text(size=6, face = "bold",color = "black", hjust=0),
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




p1 <- ggplot() +
  # Add landmarks
  geom_text(data = monterey_label, mapping = aes(x = x, y = y, label = label),
            size = 3, fontface = "bold") +
  # Add CA inset
  annotation_custom(grob = g1_inset, 
                    xmin = -122.01, 
                    xmax = -121.96,
                    ymin = 36.625) +
  #add historic kelp extent
  tidyterra::geom_spatraster(data = kelp_na, na.rm = TRUE) +
  scale_fill_gradient(low = alpha("forestgreen", alpha=0.6),
                      high = alpha("forestgreen", alpha=0.6),
                      na.value = NA,
                      labels = c("Max kelp\nextent"))+
  guides(fill = guide_legend(override.aes = list(size = 1),
                             label.theme = element_text(color = "gray"))) +
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
    data = forage_build1 %>% filter(year == 2017),
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
  labs(title = "", tag = "B")+
  theme_bw() + base_theme + theme(#axis.text.y = element_blank(),
    #legend.position = "none",
   # plot.tag.position = c(0.05, 1), 
    axis.title=element_blank()) +  
  #adjust legend position
  theme(
    legend.position = "bottom",         # Place legend at the bottom
    legend.justification = "right",     # Align the legend to the right
    legend.box = "vertical",         
    legend.margin = margin(t = -60),   # Adjust the top margin to move the legend up
    plot.margin = margin(b = 90),         # Add extra margin at the bottom of the plot
    legend.box.spacing = unit(0.001, "cm")
  )

p1


#save
#ggsave(p1, filename = file.path(figdir, "FigX_bout_locations2.png"), 
 #      width = 5, height = 4, units = "in", dpi = 600)



# Reshape the data from wide to long format
stacked_data <- sofa_build1 %>%
  dplyr::select(period, source, 2:22) %>%
  pivot_longer(cols = -c(period, source), names_to = "Prey", values_to = "Value") 

# Filter the prey types that reach or exceed 3 for mass_gain
prey_to_include <- stacked_data %>%
  filter(source == "mass_gain", Value >= 3) %>%
  distinct(Prey)

# Filter the data to include the selected prey types
filtered_data <- stacked_data %>%
  filter(Prey %in% c("urchin","mussel"))

theme2 <-  theme(axis.text=element_text(size=7, color = "black"),
                     axis.text.x = element_text(size=8,color = "black"),
                     axis.text.y = element_text(size=8,color = "black"),
                     axis.title=element_text(size=8,color = "black"),
                     legend.text=element_text(size=7,color = "black"),
                     legend.title=element_text(size=8,color = "black"),
                     plot.tag=element_text(size=10,color = "black"),
                     # Gridlines
                     panel.grid.major = element_blank(), 
                     panel.grid.minor = element_blank(),
                     panel.background = element_blank(), 
                     axis.line = element_line(colour = "black"),
                     # Legend
                     legend.key = element_rect(fill=alpha('blue', 0)),
                     legend.background = element_rect(fill=alpha('blue', 0)),
                     #facets
                     strip.text = element_text(size=6, face = "bold",color = "black", hjust=0),
                     strip.background = element_blank())


a <- ggplot(filtered_data %>% filter(period < 2020),
       aes(x = period, y = Value, color = Prey)) +
  #geom_line() +
  geom_point(alpha = 0.4) +
  #add SSW reference
  geom_vline(xintercept = 2013, linetype = "dashed", color = "black")+
  #geom_smooth(method = "loess", se = TRUE, aes(fill = Prey)) + 
  stat_smooth(geom = "line", size = 1, span = 0.7) +
  stat_smooth(method = "loess", geom = "ribbon", alpha = 0.2, aes(fill = Prey), color = NA) +
  facet_wrap(~ source, ncol = 1, scales = "free_y", 
             labeller = labeller(source = 
                                   c(dietcomp = "Proportion diet",
                                     kcal = "Energy gain (Kcal / min)",
                                     mass_gain = "Mass gain (g/min)"))) +
  labs(
    title = "",
    x = "Year",
    y = "Prey value",
    tag = "A"
  ) +
  scale_x_continuous(breaks = seq(2007, 2020, by = 2)) +  
  theme_bw() + theme2+
  theme(legend.position = "top")+
  scale_color_brewer(palette = "Dark2")+
  scale_fill_brewer(palette = "Dark2") +
  theme(plot.tag.position = c(0,0.88))
a


#g <- gridExtra::grid.arrange(a,p1, row=1)

#ggsave(g, filename = file.path(figdir, "Fig3_foraging_timeseries.png"), 
 #      width = 7, height = 9, units = "in", dpi = 600, bg = "white")


# Save the ggplot objects to grobs
a_grob <- ggplotGrob(a)
p1_grob <- ggplotGrob(p1)

# Adjust the right margin of p1_grob
p1_grob$widths[8] <- p1_grob$widths[8] + unit(1, "cm")  # Increase the margin by 1 cm

# Define the amount to nudge p1_grob down 
nudge_down <- unit(-3, "lines")

# Add a row of padding to p1_grob
p1_grob <- gtable::gtable_add_rows(p1_grob, heights = nudge_down)

# Arrange the grobs
combined_grob <- gridExtra::arrangeGrob(a_grob, p1_grob, ncol = 2, widths = c(1,1.1))

ggsave(combined_grob, filename = file.path(figdir, "Fig3_foraging_timeseries.png"), 
       width = 7, height = 6, units = "in", dpi = 600, bg = "white")


