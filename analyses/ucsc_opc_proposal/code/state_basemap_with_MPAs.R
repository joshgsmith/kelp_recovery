

# Read data
################################################################################

# Clear workspace
rm(list = ls())

# Packages
library(dplyr)
library(tidyverse)
library(patchwork)

# Directories

basedir <- "/Volumes/seaotterdb$/kelp_recovery/data"
gisdir <- file.path(basedir, "gis_data/processed")
plotdir <- here::here("analyses","ucsc_opc_proposal","figures")


# Read data
state_waters_poly <- readRDS(file.path(gisdir, "CA_state_waters_polygons.Rds"))
state_waters_line <- readRDS(file.path(gisdir, "CA_state_waters_polyline.Rds"))
mpas_orig <- readRDS(file.path(gisdir, "CA_MPA_polygons.Rds"))

# Get land
usa <- rnaturalearth::ne_states(country="United States of America", returnclass = "sf")
foreign <- rnaturalearth::ne_countries(country=c("Canada", "Mexico"), returnclass = "sf")


# Format data
################################################################################

# Reduce to MPAs of interest
mpas <- mpas_orig %>% 
  dplyr::select(name, lat_dd, long_dd, type)


# Plot data
################################################################################

# Key cities
# Monterey	36.603056	-121.893611
# Morro Bay	35.367222	-120.846667
# Pismo Beach	35.148333	-120.648056
cities <- tibble(city=c("Monterey", "Greater Farallones", "Trinidad"),
                 long_dd=c(-121.893611, -122.770086, -124.144055),
                 lat_dd=c(36.603056, 37.952678, 41.060897),
                 hjust=c(-0.2, -0.2, -0.2),
                 vjust=c(0.5, 0.5, 0.5))


# Theme
base_theme <-  theme(axis.text=element_text(size=7,colour = "black"),
                     axis.title=element_text(size=8,colour = "black"),
                     legend.text=element_text(size=7,colour = "black"),
                     legend.title=element_text(size=8,colour = "black"),
                     plot.tag=element_text(size=8,colour = "black", face = "bold"),
                     # Gridlines
                     panel.grid.major = element_blank(), 
                     panel.grid.minor = element_blank(),
                     panel.background = element_blank(), 
                     axis.line = element_line(colour = "black"),
                     # Legend
                     legend.key = element_rect(fill=alpha('blue', 0)),
                     legend.background = element_rect(fill=alpha('blue', 0)))



# Plot data
g1 <- ggplot() +
  # Plot state waters
  geom_sf(data=state_waters_line, color="grey60", lwd=0.1) +
  # Plot land
  geom_sf(data=foreign, fill="tan", color="white", lwd=0.3) +
  geom_sf(data=usa, fill="tan", color="white", lwd=0.3) +
  # Plot MPAs
  geom_sf(data=mpas, fill=ifelse(mpas$type == "SMR","indianred1","lightblue"), color="black") +
  # Plot cities
  #geom_point(data=cities, mapping=aes(x=long_dd, y=lat_dd), size=2.2) +
  #geom_text(data=cities, mapping=aes(x=long_dd, y=lat_dd, label=city, hjust=hjust, vjust=vjust),
  #          size=2.6) +
  # Labels
  labs(x="", y="") +
  #add scalebar
 # ggsn::scalebar(x.min = -122.6, x.max = -120.5, 
  #               y.min = 34.44, y.max = 37.2,
   #              #anchor=c(x=-124.7,y=41),
    #             location="bottomleft",
     #            dist = 50, dist_unit = "km",
      #           transform=TRUE, 
       #          model = "WGS84",
        #         st.dist=0.012,
         #        st.size=2,
          #       border.size=.5,
           #      height=.01
  #)+
  # Crop
  coord_sf(xlim = c(-125.6, -119), ylim = c(33.5, 41.8)) +
  # Theme
  theme_bw() + base_theme +
  theme(axis.text.y = element_text(angle = 90, hjust = 0.5),
        axis.title=element_blank(),
        legend.position = c(0.2, 0.2), 
        legend.key.size = unit(0.4, "cm")) 
g1


ggsave(g1, filename = file.path(plotdir, "CA_MPA_basemap.png"), 
       width = 8.5, height = 11, units = "in", dpi = 600)


