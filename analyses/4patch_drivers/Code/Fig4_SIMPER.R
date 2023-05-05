#Communtiy change analyses
#Joshua G. Smith; jossmith@mbayaq.org

rm(list=ls())

librarian::shelf(tidyverse, here, vegan, ggplot2, cluster, ggforce, cowplot)


################################################################################
#set directories and load data
basedir <- "/Volumes/seaotterdb$/kelp_recovery/"
figdir <- here::here("analyses","4patch_drivers","Figures")

#load raw dat
fish_raw <- read.csv(file.path(basedir, "data/subtidal_monitoring/processed/kelp_fish_counts_CC.csv")) %>%
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

upc_raw <- read.csv(file.path(basedir, "data/subtidal_monitoring/processed/kelp_upc_cov_CC.csv")) %>%
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
#calculate transect means to reduce memory for SIMPER

#drop species that were never encountered
fish_build1 <- fish_raw %>% dplyr::select(where(~ any(. != 0)))
upc_build1 <- upc_raw %>% dplyr::select(where(~ any(. != 0)))
swath_build1 <- swath_raw %>% dplyr::select(where(~ any(. != 0)))


#take mean spp counts per transect
fish_build2 <- fish_build1 %>%
  dplyr::group_by(year, MHW, baseline_region, latitude, longitude, site,
                  affiliated_mpa, mpa_class, mpa_designation) %>%
  dplyr::summarize(across(4:59, mean, na.rm = TRUE))


upc_build2 <- upc_build1 %>%
  dplyr::group_by(year, MHW, baseline_region, latitude, longitude, site,
                  affiliated_mpa, mpa_class, mpa_designation) %>%
  dplyr::summarize(across(3:41, mean, na.rm = TRUE))


swath_build2 <- swath_build1 %>%
  dplyr::group_by(year, MHW, baseline_region, latitude, longitude, site,
                  affiliated_mpa, mpa_class, mpa_designation) %>%
  dplyr::summarize(across(3:59, mean, na.rm = TRUE))

################################################################################
#process data for SIMPER

swath_build3 <- swath_build2 %>%
                        ###note::: remove these species, since these might overwhelm other changes
                        dplyr::select(!c(strongylocentrotus_purpuratus, macrocystis_pyrifera,
                                         mesocentrotus_franciscanus))
swath_groups <- swath_build3 %>% dplyr::select(1:9)
swath_dat <- swath_build3 %>% ungroup() %>% dplyr::select(10:ncol(.))


fish_groups <- fish_build2 %>% dplyr::select(1:9)
fish_dat <- fish_build2 %>% ungroup() %>% dplyr::select(10:ncol(.))

upc_groups <- upc_build2 %>% ungroup() %>% dplyr::select(1:9)
upc_dat <- upc_build2 %>% ungroup() %>% dplyr::select(10:ncol(.))

################################################################################

#run SIMPER across MHW periods

sim_swath <- with(swath_groups, simper(swath_dat, MHW), 
                  ordered=TRUE)

sim_fish <- with(fish_groups, simper(fish_dat, MHW), 
                  ordered=TRUE)

sim_upc <- with(upc_groups, simper(upc_dat, MHW, parallel = TRUE), 
                 ordered=TRUE)

################################################################################
#collect output

swath_table <- as.data.frame(summary(sim_swath)$before_after)%>%
  mutate(method="swath",
         contrib = cumsum-lag(cumsum, default=0),
         perc_change = (avb-ava)/ava,
         sign = ifelse(perc_change > 0, "positive","negative"))

fish_table <- as.data.frame(summary(sim_fish)$before_after)%>%
  mutate(method="fish",
         contrib = cumsum-lag(cumsum, default=0),
         perc_change = (avb-ava)/ava,
         sign = ifelse(perc_change > 0, "positive","negative"))

upc_table <- as.data.frame(summary(sim_upc)$before_after)%>%
  mutate(method="UPC",
         contrib = cumsum-lag(cumsum, default=0),
         perc_change = (avb-ava)/ava,
         sign = ifelse(perc_change > 0, "positive","negative"))



simper_table <- rbind(swath_table, 
                          fish_table,
                          upc_table)%>%
  tibble::rownames_to_column(var = "species") %>%
  filter(cumsum <0.80)%>%
  dplyr::select(method, species, avg_before = avb, avg_after=ava, cumsum, contrib, perc_change, sign) %>%
  #filter to individual contributions > 5%
  filter(contrib > 0.05)


################################################################################
#select top contributors from raw data

swath_plot <- swath_raw %>% dplyr::select(year,MHW, site, zone, transect,
                                          macrocystis_pyrifera, strongylocentrotus_purpuratus, mesocentrotus_franciscanus,
                                          #balanus_nubilus,
                                          pterygophora_californica,
                                          laminaria_setchellii) %>%
                pivot_longer(cols = c(macrocystis_pyrifera, strongylocentrotus_purpuratus, mesocentrotus_franciscanus,
                                      #balanus_nubilus,
                                      pterygophora_californica,
                                      laminaria_setchellii), names_to = "species",
                             values_to = "counts") %>%
                mutate(MHW = str_to_title(factor(MHW, levels=c('before','during','after'))),
                       species = ifelse(species == "balanus_nubilus","Balanus nubilus",ifelse(
                         species == "laminaria_setchellii","Laminaria setchellii",ifelse(
                           species == "macrocystis_pyrifera","Macrocystis pyrifera",ifelse(
                             species == "mesocentrotus_franciscanus","Mesocentrotus franciscanus",ifelse(
                               species == "pterygophora_californica","Pterygophora californica",ifelse(
                                 species == "strongylocentrotus_purpuratus","Strongylocentrotus purpuratus",species
                               )
                             )
                           )
                         )
                       )),
                       `Organism type` = factor(ifelse(species == "Macrocystis pyrifera" | species == "Pterygophora californica" | species == "Laminaria setchellii","Stipitate algae",ifelse(
                         species == "Strongylocentrotus purpuratus","Purple sea urchin","Red sea urchin"
                       ))))
swath_plot$MHW <- factor(swath_plot$MHW, levels = c("Before", "During", "After"))
swath_plot$species <- factor(swath_plot$species, levels = c("Strongylocentrotus purpuratus","Mesocentrotus franciscanus","Macrocystis pyrifera","Pterygophora californica","Laminaria setchellii"))

fish_plot <- fish_raw %>% dplyr::select(year,MHW, site, zone,level, transect,
                                        sebastes_mystinus,
                                        oxyjulis_californica,
                                        sebastes_serranoides_flavidus) %>%
  pivot_longer(cols = c(sebastes_mystinus,
                        oxyjulis_californica,
                        sebastes_serranoides_flavidus), names_to = "species",
               values_to = "counts") %>%
  mutate(MHW = str_to_title(factor(MHW, levels=c('before','during','after'))),
         species = ifelse(species == "sebastes_mystinus","Sebastes mystinus",ifelse(
           species == "oxyjulis_californica","Oxyjulis californica","Sebastes serranoides"
         )),
         `Organism type` = "Fish")
swath_plot$MHW <- factor(swath_plot$MHW, levels = c("Before", "During", "After"))


upc_plot <- upc_raw %>% dplyr::select(year,MHW, site, zone, transect,
                                      coralline_algae_crustose, 
                                      #coralline_algae_erect_articulated,
                                      red_algae_branching_flat_blade, red_algae_encrusting,
                                      red_algae_leaf_like) %>%
                   pivot_longer(cols = c(coralline_algae_crustose, 
                                         #coralline_algae_erect_articulated,
                                         red_algae_branching_flat_blade, red_algae_encrusting,
                                         red_algae_leaf_like), names_to = "species",
                            values_to = "counts") %>%
                  #rename for plotting
                  mutate(species = ifelse(species == "coralline_algae_crustose","Crustose coralline algae",ifelse(
                    species == "coralline_algae_erect_articulated","Articulated coralline algae",ifelse(
                      species == "red_algae_branching_flat_blade","Flat blade red macroalgae",ifelse(
                        species == "red_algae_encrusting","Encrusting red algae",ifelse(
                          species == "red_algae_leaf_like","Leaf-like red macroalgae",species
                        )
                      )
                    )
                  ))) %>%
                  mutate(MHW = str_to_title(factor(MHW, levels=c('before','during','after'))),
                         `Organism type` = ifelse(species == "Flat blade red macroalgae"| species == "Leaf-like red macroalgae","Macroalgae","Encrusting"))
upc_plot$MHW <- factor(upc_plot$MHW, levels = c("Before", "During", "After"))
upc_plot$species <- factor(upc_plot$species, levels = c("Flat blade red macroalgae","Leaf-like red macroalgae", "Crustose coralline algae","Encrusting red algae"))
upc_plot$`Organism type` <- factor(upc_plot$`Organism type`, levels = c("Macroalgae","Encrusting"))

################################################################################

################################################################################
#plot

# Schematic theme


my_theme <-  theme(axis.text=element_text(size=6),
                   axis.text.y = element_text(angle = 90, hjust = 0.5),
                   axis.title=element_text(size=8),
                   plot.tag=element_text(size=8, face="bold"),
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



p1 <- ggplot(upc_plot, aes(x = MHW, y = counts, fill = `Organism type`)) +
  geom_boxplot() +
  facet_wrap(~ species, scales = "fixed", nrow=1) +
  scale_fill_manual(name = "Organism type",
    values= c("Macroalgae"='#3CB371',"Encrusting" = 'pink', "Purple sea urchin" = "purple", "Red sea urchin" = "indianred","Stipitate algae" = "brown",
              drop=FALSE))+
  theme_bw() +
  my_theme+
  labs(tag = "A", fill = "Organism type")+
  ylab("Percent cover")+
  xlab("Heatwave period")
p1


p2 <- ggplot(swath_plot, aes(x = MHW, y = log(counts))) +
  geom_boxplot( aes(fill = `Organism type`)) +
  facet_wrap(~ species, scales = "free_y", nrow=1) +
  #scale_fill_manual(values= c('#3CB371','pink'))+
  theme_bw() +
  my_theme+
  scale_fill_manual(name = "Organism type",
                    values= c("Macroalgae"='#3CB371',"Encrusting" = 'pink', "Purple sea urchin" = "purple", "Red sea urchin" = "indianred","Stipitate algae" = "brown",
                              drop=FALSE))+
  labs(tag = "B", fill = "")+
  ylab("Log density \n(no. individuals per 60 m²)")+
  xlab("Heatwave period")
p2


p3 <- ggplot(fish_plot, aes(x = MHW, y = log(counts), fill = `Organism type`)) +
  geom_boxplot( aes(fill = `Organism type`)) +
  facet_wrap(~ species, scales = "free_y", nrow=1) +
  theme_bw() +
  my_theme+
  scale_fill_manual(name = "Organism type",
                    values= c("Macroalgae"='#3CB371',"Encrusting" = 'pink', "Purple sea urchin" = "purple", "Red sea urchin" = "indianred","Stipitate algae" = "brown",
                              "Fish" = "lightblue",
                              drop=FALSE))+
  labs(tag = "C", fill = "")+
  ylab("Log density \n(no. individuals per 60 m²)")+
  xlab("Heatwave period")
p3


combined_plot <- ggpubr::ggarrange(p1,p2,p3,ncol=1, common.legend = TRUE, legend = "right")
combined_plot




library(patchwork)

combined <- p1 + p2 + p3 + plot_layout(
  ncol=1,
  guides = "collect") + theme(legend.position = "right",
                              legend.spacing = unit(-4, "cm"))
combined


ggsave(combined, filename=file.path(figdir, "Fig4_boxplot.png"), bg = "white",
       width=7, height=8, units="in", dpi=600) 






