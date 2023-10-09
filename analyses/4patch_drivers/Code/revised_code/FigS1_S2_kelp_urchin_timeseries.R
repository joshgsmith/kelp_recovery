#Community data processing
#Joshua G. Smith
#May 1, 2023

rm(list=ls())
librarian::shelf(tidyverse, here, vegan)


################################################################################
#set directories and load data
basedir <- "/Volumes/seaotterdb$/kelp_recovery/"
figdir <- here::here("analyses","4patch_drivers","Figures")

swath_raw <- read.csv(file.path(basedir, "data/subtidal_monitoring/processed/kelp_swath_counts_CC.csv")) 

################################################################################
#select and reshape data

swath_sub1 <- swath_raw %>% dplyr::select(1:11, 'macrocystis_pyrifera', 'strongylocentrotus_purpuratus') %>%
                            pivot_longer(cols = c('macrocystis_pyrifera', 'strongylocentrotus_purpuratus'),
                                         names_to = "species", values_to="counts") %>%
  mutate(outbreak_period = ifelse(year <2014, "before","after"))%>%
  dplyr::select(year, outbreak_period, everything())%>%
  data.frame() 

#calculate mean kelp density pre-2013

kelp_mean <- swath_sub1 %>% filter(species == "macrocystis_pyrifera") %>% group_by(year, site, outbreak_period)%>%
                summarise(kelp_mean = mean(counts, na.rm=TRUE),
                          one_sd = sd(counts, na.rm=TRUE)) 

                

#join
swath_sub <- left_join(swath_sub1, kelp_mean, by=c("year", "site", "outbreak_period"))

################################################################################
#plot distinct sites

site_year <- swath_raw %>% dplyr::select(year, site) %>% distinct() %>%
  mutate(
  site_order = case_when(
    site == "CANNERY_UC" ~ 1,
    site == "CANNERY_DC" ~ 2,
    site == "MACABEE_UC" ~ 3,
    site == "MACABEE_DC" ~ 4,
    site == "HOPKINS_UC" ~ 5,
    site == "HOPKINS_DC" ~ 6,
    site == "LOVERS_UC" ~ 7,
    site == "LOVERS_DC" ~ 8,
    site == "SIREN" ~ 9,
    site == "OTTER_PT_UC" ~ 10,
    site == "OTTER_PT_DC" ~ 11,
    site == "LONE_TREE" ~ 12,
    site == "PESCADERO_UC" ~ 13,
    site == "PESCADERO_DC" ~ 14,
    site == "STILLWATER_UC" ~ 15,
    site == "STILLWATER_DC" ~ 16,
    site == "BUTTERFLY_UC" ~ 17,
    site == "BUTTERFLY_DC" ~ 18,
    site == "MONASTERY_UC" ~ 19,
    site == "MONASTERY_DC" ~ 20,
    site == "BLUEFISH_UC" ~ 21,
    site == "BLUEFISH_DC" ~ 22,
    site == "WESTON_UC" ~ 23,
    site == "WESTON_DC" ~ 24,
    TRUE ~ NA)
)

# Theme
my_theme <- theme(
  axis.text = element_text(size = 6),
  axis.text.x = element_text(angle = 45, hjust = 1),
  axis.title = element_text(size = 8),
  plot.tag = element_blank(),
  plot.title = element_text(size = 7, face = "bold"),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  panel.background = element_blank(),
  axis.line = element_line(colour = "black"),
  legend.key = element_blank(),
  legend.background = element_rect(fill = alpha('blue', 0)),
  legend.key.height = unit(1, "lines"),
  legend.text = element_text(size = 6),
  legend.title = element_text(size = 7),
  strip.background = element_blank(),
  strip.text = element_text(size = 6, face = "bold"),
  # Adjust panel spacing
  panel.spacing.y = unit(0.5, "lines")
)

g <- ggplot(site_year %>% mutate(site = str_to_title(gsub("_", " ", site)),
                            site = str_replace(site, "Dc", "DC"),
                            site = str_replace(site, "Uc", "UC")),
       aes(x = year, y = reorder(site, site_order))) +
  geom_tile(fill = "#8B5FBF", colour = "white", alpha=0.7) +
  labs(x = "Year", y = "Site") +
  scale_x_continuous(breaks = unique(site_year$year), labels = unique(site_year$year)) +
  theme_bw()+
  my_theme
g

#ggsave(g, filename=file.path(figdir, "FigS1_site_frequency.png"), 
 #      width=6, height=6, units="in", dpi=600)

################################################################################
#plot

# Theme
my_theme <-  theme(axis.text=element_text(size=6),
                   axis.text.y = element_text(angle = 90, hjust = 0.5),
                   axis.title=element_text(size=8),
                   plot.tag=element_blank(), #element_text(size=8),
                   plot.title =element_text(size=7, face="bold"),
                   # Gridlines 
                   panel.grid.major = element_blank(), 
                   panel.grid.minor = element_blank(),
                   panel.background = element_blank(), 
                   axis.line = element_line(colour = "black"),
                   # Legend
                   legend.key = element_blank(),
                   legend.background = element_rect(fill=alpha('blue', 0)),
                   legend.key.height = unit(1, "lines"), 
                   legend.text = element_text(size = 6),
                   legend.title = element_text(size = 7),
                   #legend.spacing.y = unit(0.75, "cm"),
                   #facets
                   strip.background = element_blank(),
                   strip.text = element_text(size = 6 ,face="bold"),
                   )


# create a new column for the site factor with levels ordered by kelp_mean


g <- ggplot(swath_sub %>%
        mutate(species = ifelse(species == "macrocystis_pyrifera","Giant kelp \n(M. pyrifera)","Purple sea urchins \n(S. purpuratus)"),
               site = gsub("_", " ", site))
       , aes(x = year, y = pmax(log(counts), 0.001), fill = species, color=species, group=species)) +
  geom_point(aes(color = species), alpha = 0.2, size=0.5) +
  geom_smooth(method = "auto", se = TRUE, size = 0.5, aes(group = species, color = species), alpha = 0.3, inherit.aes = TRUE) +
  #SSW
  geom_vline(xintercept = 2013, linetype = "dotted", size=0.3)+
  annotate(geom="text", label="SSW", x=2010, y=10 , size=2) +
  annotate("segment", x = 2011.8, y = 9.8, xend = 2012.7, yend = 8.8,
           arrow = arrow(type = "closed", length = unit(0.02, "npc")))+
  #define MHW period
  # Heatwave
  annotate(geom="rect", xmin=2014, xmax=2016, ymin=-Inf, ymax=Inf, fill="indianred1", alpha=0.2) +
  annotate(geom="text", label="MHW", x=2018.5, y=10 , size=2) +
  annotate("segment", x = 2016.5, y = 10, xend = 2015, yend = 10,
           arrow = arrow(type = "closed", length = unit(0.02, "npc")))+
  #
  facet_wrap(~ site, scales = "free") +
  scale_color_manual(values = c("forestgreen", "purple")) +
  scale_fill_manual(values = c("forestgreen", "purple")) +
  ylim(-1, 10) +
  labs(fill = "Species", color = "Species")+
  ylab("Log density (No. individuals per 60 m²)")+
  xlab("Year")+
  my_theme + theme(legend.position = "top")

g

# Export figure
ggsave(g, filename=file.path(figdir, "FigS1_urchin_kelp_timeseries.png"), 
      width=6.5, height=7, units="in", dpi=600)

################################################################################
#plot region
# Theme
my_theme <-  theme(axis.text=element_text(size=8),
                   axis.text.y = element_text(angle = 90, hjust = 0.5),
                   axis.title=element_text(size=8),
                   plot.tag=element_blank(), #element_text(size=8),
                   plot.title =element_text(size=9, face="bold"),
                   # Gridlines 
                   panel.grid.major = element_blank(), 
                   panel.grid.minor = element_blank(),
                   panel.background = element_blank(), 
                   axis.line = element_line(colour = "black"),
                   # Legend
                   legend.key = element_blank(),
                   legend.background = element_rect(fill=alpha('white', 0)),
                   legend.key.height = unit(1, "lines"), 
                   legend.text = element_text(size = 8),
                   legend.title = element_text(size = 9),
                   #legend.spacing.y = unit(0.75, "cm"),
                   #facets
                   strip.background = element_blank(),
                   strip.text = element_text(size = 8 ,face="bold"),
)

swath_sub_site <- swath_sub %>% group_by(year, site, species) %>% summarize(count_mean = mean(counts)) %>%
  mutate(log_den = ifelse(species == "strongylocentrotus_purpuratus", log(count_mean),count_mean),
         sqrt_den = ifelse(species == "strongylocentrotus_purpuratus", sqrt(count_mean),count_mean)) %>%
  dplyr::select(year, site, species, sqrt_den)%>%
        pivot_wider(values_from = sqrt_den, names_from = species) %>%
  mutate(sqrt_5 = strongylocentrotus_purpuratus*5)


#trick ggplot to think it is monotone
f <- Vectorize(function(x) {
  if (x < 1) return(x/1e10)
  (x/5)^2 #reverse operation of sqrt transformation
})

g1 <- ggplot(swath_sub_site %>%
              mutate(
                     site = gsub("_", " ", site)), aes(x=year)) +
  geom_point(aes(y = macrocystis_pyrifera, color = "Kelp \n(M. pyrifera)"), alpha = 0.2, size=0.5) +
  geom_smooth(method = "auto", se = TRUE, size = 0.5, aes(y=macrocystis_pyrifera,color = "Kelp \n(M. pyrifera)",
                                                          fill = "Kelp \n(M. pyrifera)"
                                                          ), alpha = 0.3, inherit.aes = TRUE) +
  geom_point(aes(y = sqrt_5, color = "Sea urchins \n(S. purpuratus)"), alpha = 0.2, size=0.5) +
  geom_smooth(method = "auto", se = TRUE, size = 0.5, aes(y=sqrt_5, color = "Sea urchins \n(S. purpuratus)",
                                                          fill = "Sea urchins \n(S. purpuratus)"), alpha = 0.3, inherit.aes = TRUE) +
  #SSW
  geom_vline(xintercept = 2013, linetype = "dotted", size=0.3)+
  annotate(geom="text", label="SSW", x=2011, y=200 , size=2) +
  annotate("segment", x = 2011.8, y = 198, xend = 2012.7, yend = 188,
           arrow = arrow(type = "closed", length = unit(0.02, "npc")))+
  # Heatwave
  annotate(geom="rect", xmin=2014, xmax=2016, ymin=-Inf, ymax=Inf, fill="indianred1", alpha=0.2) +
  annotate(geom="text", label="MHW", x=2017.5, y=200 , size=2) +
  annotate("segment", x = 2016.5, y = 200, xend = 2015, yend = 200,
           arrow = arrow(type = "closed", length = unit(0.02, "npc")))+
  scale_color_manual(values = c("Kelp \n(M. pyrifera)" = "forestgreen", "Sea urchins \n(S. purpuratus)" = "purple")) +
  scale_fill_manual(values = c("Kelp \n(M. pyrifera)" = "forestgreen", "Sea urchins \n(S. purpuratus)" = "purple")) +
  scale_y_continuous(name="Macrocystis pyrifera \n(no. stipe per 60m²)", sec.axis = sec_axis(~f(.), name = "Strongylocentrotus purpuratus \n(no. per 60m²)",
                     breaks = c(2, 50, 200, 400, 800, 1600, 3200)),
                     limits = c(-3,200) #set limit of y1
                     ) +
  
  scale_fill_manual(values = c("forestgreen", "purple")) +
  #ylim(0, 300) +
  labs(fill = "Species", color = "Species")+
  guides(color = guide_legend(title = NULL),
         fill = guide_legend(title = NULL))+ 
  ylab("Log density (No. individuals per 60 m²)")+
  xlab("Year")+
  my_theme + theme(legend.position = "top",
                  axis.text.y = element_text(angle = 0, hjust = 1))

g1




# Italicize the legend labels using expression
g1 <- g1 +
  scale_color_manual(
    values = c("Kelp \n(M. pyrifera)" = "forestgreen", "Sea urchins \n(S. purpuratus)" = "purple"),
    labels = c(expression("Kelp" ~ italic("(M. pyrifera)")), expression("Sea urchins"~italic("(S. purpuratus)")))
  ) +
  scale_fill_manual(
    values = c("Kelp \n(M. pyrifera)" = "forestgreen", "Sea urchins \n(S. purpuratus)" = "purple"),
    labels = c(expression("Kelp" ~ italic("(M. pyrifera)")), expression("Sea urchins"~italic("(S. purpuratus)")))
  )
g1


# Italicize the axis labels using expression with manual line breaks
g1 <- g1 +
  scale_y_continuous(
    name = expression(italic("M. pyrifera") ~ "(stipes per 60m²)"),
    sec.axis = sec_axis(~f(.), name = expression(italic("S. purpuratus") ~ "(no. per 60m²)"),
                        breaks = c(2, 50, 200, 400, 800, 1600, 3200)),
    limits = c(-3, 200) # Set limit of y1
  ) +
  scale_x_continuous(
    name = expression("Year")
  )



