#Community data processing
#Joshua G. Smith
#May 1, 2023

rm(list=ls())
librarian::shelf(tidyverse, here, vegan)

library(dplyr)

################################################################################
#set directories and load data
basedir <- "/Volumes/seaotterdb$/kelp_recovery/"
figdir <- here::here("analyses","4patch_drivers","Figures")

swath_raw <- read.csv(file.path(basedir, "data/subtidal_monitoring/processed/kelp_swath_counts_CC.csv")) %>%
  #select sites in Carmel and Monterey Bay only
  dplyr::filter(latitude >= 36.46575 & latitude <= 36.64045) %>%
  #drop sites with insufficient data
  dplyr::filter(!(site == "ASILOMAR_DC" |
                    site == "ASILOMAR_UC" |
                    site == "CHINA_ROCK" |
                    site == "CYPRESS_PT_DC" |
                    site == "CYPRESS_PT_UC" |
                    site == "PINNACLES_IN" |
                    site == "PINNACLES_OUT" |
                    site == "PT_JOE" |
                    site == "SPANISH_BAY_DC" |
                    site == "SPANISH_BAY_UC" |
                    site == "BIRD_ROCK"))

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
                summarise(kelp_mean = mean(counts))
                


#join
swath_sub <- left_join(swath_sub1, kelp_mean, by=c("year", "site", "outbreak_period"))


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



swath_sub$counts[swath_sub$species == "macrocystis_pyrifera"] <- 
  swath_sub$counts[swath_sub$species == "macrocystis_pyrifera"] + 0.1

swath_sub$counts[swath_sub$species == "strongylocentrotus_purpuratus"] <- 
  swath_sub$counts[swath_sub$species == "strongylocentrotus_purpuratus"] + 0.1


# create a new column for the site factor with levels ordered by kelp_mean
swath_sub$site_ordered <- reorder(swath_sub$site, -swath_sub$kelp_mean)

g <- ggplot(swath_sub %>%
        mutate(species = ifelse(species == "macrocystis_pyrifera","Giant kelp \n(M. pyrifera)","Purple sea urchins \n(S. purpuratus)"),
               site = gsub("_", " ", site))
       , aes(x = year, y = log(counts), fill = species, color=species, group=species)) +
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
  ylim(-0.5, 10) +
  labs(fill = "Species", color = "Species")+
  ylab("Log density (No. individuals per 60 m²)")+
  xlab("Year")+
  my_theme + theme(legend.position = "top")

# Export figure
ggsave(g, filename=file.path(figdir, "FigS1_urchin_kelp_timeseries.png"), 
       width=6.5, height=7, units="in", dpi=600)











