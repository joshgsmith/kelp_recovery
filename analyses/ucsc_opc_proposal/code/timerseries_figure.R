

#Joshua G. Smith; jossmith@mbayaq.org

rm(list=ls())

######
######
#required packages
librarian::shelf(tidyverse, sf, raster, shiny, tmap, terra, tidyterra, RColorBrewer)

#set directories 
basedir <- "/Volumes/seaotterdb$/kelp_recovery/data"
figdir <- here::here("analyses","ucsc_opc_proposal","figures")

#read state
ca_counties_orig <- st_read(file.path(basedir, "gis_data/raw/ca_county_boundaries/s7vc7n.shp")) 

#read landsat raw
landsat_orig <- st_read(file.path(basedir, "kelp_landsat/processed/monterey_peninsula/landsat_mpen_1984_2023_points_withNAs.shp"))

#read subtidal data
swath_orig <- read_csv(file.path(basedir,"subtidal_monitoring/processed/kelp_swath_raw_CC.csv"))

#load envr data
mod_predict <- readRDS(file.path(basedir, "/environmental_data/predictors_at_pisco_sites_v2.Rds"))

# Get land
usa <- rnaturalearth::ne_states(country="United States of America", returnclass = "sf")
foreign <- rnaturalearth::ne_countries(country=c("Canada", "Mexico"), returnclass = "sf")


################################################################################
# summarize landsat data to determine baseline avg and standard deviations

# Group by year and quarter and summarize the biomass values
# We want the total annual bioamss for the region
summarized_data <- landsat_orig %>%
  #filter to year 1990 and beyond for Q3 only
  filter(year >= 1985 & quarter == 3)%>%
  group_by(year) %>%
  summarize(total_biomass = sum(biomass, na.rm = TRUE))

#determine baseline average kelp biomass
baseline_average <- summarized_data %>%
  filter(year < 2014) %>%
  summarize(baseline_avg = mean(total_biomass, na.rm=TRUE),
            baseline_sd = sd(total_biomass, na.rm=TRUE)) 

#inspect
#View(baseline_average)

# calculate departures in standard deviation from baseline avg
landsat_data <- summarized_data %>%
  st_join(baseline_average) %>%
  mutate(deviation = (total_biomass - baseline_avg) / baseline_sd)

################################################################################
# Determine baseline mean and deviation for urchins and pycno

# Group by year and take the mean counts per transect
summarized_data <- swath_orig %>%
  # Drop sites
  dplyr::filter(!(affiliated_mpa == "Natural Bridges SMR" | affiliated_mpa == "Vandenberg SMR" | affiliated_mpa == "Cambria SMCA")) %>%
  dplyr::select(1:11, strongylocentrotus_purpuratus, pycnopodia_helianthoides) %>%
  # Make longer
  pivot_longer(cols = c(strongylocentrotus_purpuratus, pycnopodia_helianthoides), 
               names_to = "species", 
               values_to = "counts") %>%
  filter(year >= 2000) %>%
  group_by(year, species) %>%
  summarize(mean_counts = mean(counts, na.rm = TRUE))

# Determine baseline average and standard deviation for each species
baseline_average <- summarized_data %>%
  group_by(species) %>%
  filter(year < 2014) %>%
  summarize(baseline_avg = mean(mean_counts, na.rm = TRUE),
            baseline_sd = sd(mean_counts, na.rm = TRUE)) 

# Calculate departures in standard deviation from baseline avg
swath_data <- summarized_data %>%
  left_join(baseline_average) %>%
  mutate(deviation = (mean_counts - baseline_avg) / baseline_sd)

################################################################################
# process SST dat

sst_data <- mod_predict %>% group_by(year) %>% summarize(mean_anom = mean(sst_month_anom))



################################################################################
# plot

# Theme
base_theme <-  theme(axis.text=element_text(size=12, color = "black"),
                     axis.title=element_text(size=6,color = "black"),
                     plot.tag=element_text(size=7,color = "black"),
                     plot.title=element_text(size=6,color = "black", face = "bold"),
                     # Gridlines
                     panel.grid.major = element_blank(), 
                     panel.grid.minor = element_blank(),
                     panel.background = element_blank(), 
                     axis.line = element_line(colour = "black"),
                     # Legend
                     legend.key.size = unit(0.3, "cm"), 
                     legend.key = element_rect(fill=alpha('blue', 0)),
                     legend.spacing.y = unit(0.1, "cm"),  
                     legend.text=element_text(size=4,color = "black"),
                     legend.title=element_text(size=4,color = "black"),
                     #legend.key.height = unit(0.1, "cm"),
                     legend.background = element_rect(fill=alpha('blue', 0)),
                     #facets
                     strip.text = element_text(size=6, face = "bold",color = "black", hjust=0),
                     strip.background = element_blank())


# Calculate the scaling factor to fit the positive values between 0 and 2.5
positive_scaling_factor <- 2.5 / max(swath_data$deviation[swath_data$deviation >= 0])

# Center and scale the positive deviation values
swath_data$transformed_deviation <- ifelse(swath_data$species == "strongylocentrotus_purpuratus" & swath_data$deviation >= 0,
                                           swath_data$deviation * positive_scaling_factor,
                                           swath_data$deviation)

# Create the ggplot object using the final_data dataset
p <- ggplot() +
  # Add swath data and create stacked barplots with transformed y values
  geom_bar(aes(x = year, 
               y = transformed_deviation,
               fill = species), 
           stat = "identity", position = "stack", width = 0.8, data = swath_data) +
  # Set colors for the species with conditional alpha
  scale_fill_manual(
    values = c("strongylocentrotus_purpuratus" = "#663398", "pycnopodia_helianthoides" = alpha("pink", 0.6))) +
  geom_hline(yintercept = 0, linetype = "dotted", color = "black") +
  labs(x = "Year", y = "Transformed Deviation") +  # Update the y-axis label
  scale_x_continuous(breaks = seq(1985, 2021, by = 5)) +
  annotate("text", x = 2021, y = Inf, label = "b", color = "black", hjust = 1, vjust = 1) +
  # Add landsat
  geom_point(data = subset(landsat_data, year <= 2021), aes(x = year, y = pmin(deviation, 3)), color = "brown") +
  geom_line(data = subset(landsat_data, year <= 2021), aes(x = year, y = pmin(deviation, 3)), color = "brown") +
  # Add sst
  geom_point(data = subset(sst_data, year <= 2021), aes(x = year, y = mean_anom, color = "orange")) +
  geom_line(data = subset(sst_data, year <= 2021), aes(x = year, y = mean_anom), color = "orange") +
  scale_y_continuous(limits = c(-2.5, 3)) +
  theme_bw() + base_theme +
  theme(legend.position = "none")  # Hide the legend

# Display the plot
print(p)

ggsave(p, filename = file.path(figdir, "Fig1_timeseries.png"), 
     width = 8, height = 2.5, units = "in", dpi = 600)



