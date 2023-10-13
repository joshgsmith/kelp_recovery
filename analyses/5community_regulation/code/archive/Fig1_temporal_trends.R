
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



g <- ggplot(combined_data, aes(x = year, y = response)) +
  geom_point(aes(color = category), alpha = 0.4) +
  #stars
  stat_smooth(data = combined_data %>% filter(category == "P. ochraceus"), geom = "line", size = 1, span = 0.6, aes(color = category)) +
  stat_smooth(data = combined_data %>% filter(category == "P. ochraceus"), method = "loess", geom = "ribbon", alpha = 0.2, aes(fill = category), color = NA, span = 0.6) +
  #mussels
  stat_smooth(data = combined_data %>% filter(category == "Mussels"), geom = "line", size = 1, span = 0.6, aes(color = category)) +
  stat_smooth(data = combined_data %>% filter(category == "Mussels"), method = "loess", geom = "ribbon", alpha = 0.2, aes(fill = category), color = NA, span = 0.6) +
  #otters
  stat_smooth(data = combined_data %>% filter(category == "Total independent otters"), geom = "line", size = 1, span = 0.4, aes(color = category)) +
  stat_smooth(data = combined_data %>% filter(category == "Total independent otters"), method = "loess", geom = "ribbon", alpha = 0.2, aes(fill = category), color = NA, span = 0.4) +
  facet_wrap(~category, scales = "free_y", ncol = 1) +
  geom_vline(xintercept = 2013, linetype = "dotted", size = 0.6) +
  labs(x = "Year", y = "Total Count") +
  scale_color_brewer(palette = "Dark2", guide = guide_legend()) +  # Use Dark2 palette
  ggtitle("") +
  theme_minimal() +
  theme(legend.title = element_blank()) + theme_bw() + base_theme+
  scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.05)), oob = scales::squish)

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
  labs(x = "Year", y = "Density (no. per m²)", title = "P. ochraceus") 

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
  labs(x = "", y = "Percent cover", title = "Mussels") +
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
  labs(x = "", y = "Number of independents",title = "Sea otters") +
  theme(axis.text.x = element_blank())

p3


p <- gridExtra::grid.arrange(p3,p2,p1, ncol=1)


#save
ggsave(p, filename = file.path(figdir, "Fig1_timeseries.png"), 
       width =5, height = 6, units = "in", dpi = 600)


