#Community data processing
#Joshua G. Smith
#May 1, 2023

rm(list=ls())
librarian::shelf(tidyverse, here, vegan, reshape2)


################################################################################
#set directories and load data
basedir <- "/Volumes/seaotterdb$/kelp_recovery/"
figdir <- here::here("analyses","4patch_drivers","Figures")
outdir <- here::here("analyses","4patch_drivers","Output")

load(file.path(outdir,"multivariate_data.Rdata"))

swath_raw <- read.csv(file.path(basedir, "data/subtidal_monitoring/processed/kelp_swath_counts_CC.csv")) 

################################################################################
#generate distance matrix

set.seed(1985)

#define group vars
stan_group_vars <- stan_dat %>% dplyr::select(year, site) %>%
  #define urchin outbreak period
  mutate(period = ifelse(year > 2013, "after","before"),
         identifier = paste(year, site, period)) %>%
  dplyr::select(identifier)

#define data for matrix
stan_ord_dat <- stan_dat %>% dplyr::select(10:ncol(.))

#create matrix
stan_distmat <- as.matrix(vegdist(stan_ord_dat, method = "bray", na.rm = T))

#join group vars
dist_dat <- cbind(stan_group_vars, stan_distmat)

#create header names to match square matrix
colnames(dist_dat)[2:ncol(dist_dat)] <- dist_dat[,1]

#convert to three column format
dist_long <- setNames(melt(dist_dat), c('site_year_period1', 'site_year_period2', 'dissim')) 

dist_long1 <- dist_long %>% 
  #bust apart the grouping var
  mutate(year_1 = word(site_year_period1,1),
         year_2 = word(site_year_period2,1),
         period_1 = word(site_year_period1,-1),
         period_2 = word(site_year_period2,-1),
         site_1 = word(site_year_period1,-2),
         site_2 = word(site_year_period2,-2)) %>%
  #filter site-by-site comparisons only
  filter(site_1==site_2,
         period_1 == "after",
         period_2 == "before"
  ) %>%
  mutate(sim = 1-dissim)%>%
  #drop duplicates
  distinct()%>%
  #caclulate mean sim per site
  group_by(site_1)%>%
  summarize(mean_sim = mean(sim)) %>%
  rename(site = site_1)



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
  summarise(kelp_mean = mean(counts, na.rm=TRUE))



#join
swath_sub <- left_join(swath_sub1, kelp_mean, by=c("year", "site", "outbreak_period"))


################################################################################
#calculate % change of kelp and urchins at the transect level

df_mean_kelp <- swath_sub %>%
  group_by(site, species, outbreak_period) %>% 
  dplyr::summarise(mean_count = mean(counts, na.rm=TRUE),
                   se_count = sd(counts, na.rm = TRUE)/sqrt(sum(!is.na(counts))))%>%
  pivot_wider(names_from = outbreak_period, values_from = c(mean_count, se_count)) %>%
  mutate_all(~ if_else(is.na(.), 0, .))


################################################################################
#join with sim 

plot_dat <- left_join(df_mean_kelp, dist_long1, by="site")

################################################################################
#determine slopes

slope_dat <- plot_dat %>% dplyr::select(1:4, mean_sim) %>% 
  pivot_wider(names_from = species, values_from = c(mean_count_before, mean_count_after)) %>%
  #calculate slope
  mutate(slope = (mean_count_after_macrocystis_pyrifera - mean_count_before_macrocystis_pyrifera) /
                 (mean_count_after_strongylocentrotus_purpuratus - mean_count_before_strongylocentrotus_purpuratus),
         perc_change_kelp = (mean_count_after_macrocystis_pyrifera - mean_count_before_macrocystis_pyrifera)/mean_count_before_macrocystis_pyrifera,
         perc_change_urch =  (mean_count_after_strongylocentrotus_purpuratus - mean_count_before_strongylocentrotus_purpuratus)/ mean_count_before_strongylocentrotus_purpuratus)
           


# Fit a linear regression model
model <- lm(mean_sim ~ slope, data = slope_dat)
model <- lm(mean_sim ~ perc_change_kelp, data = slope_dat)
model <- lm(mean_sim ~ perc_change_urch, data = slope_dat)

summary(model)


################################################################################
#plot

# Theme
my_theme <-  theme(axis.text=element_text(size=8, color = "black"),
                   axis.text.y = element_text(angle = 90, hjust = 0.5, color ="black"),
                   axis.title=element_text(size=10, color = "black"),
                   plot.tag= element_text(size=8, color = "black"),
                   plot.title =element_text(size=10, face="italic", color = "black"),
                   # Gridlines 
                   panel.grid.major = element_blank(), 
                   panel.grid.minor = element_blank(),
                   panel.background = element_blank(), 
                   axis.line = element_line(colour = "black"),
                   # Legend
                   legend.key = element_blank(),
                   legend.background = element_rect(fill=alpha('blue', 0)),
                   legend.key.height = unit(1, "lines"), 
                   legend.text = element_text(size = 6, color = "black"),
                   legend.title = element_text(size = 7, color = "black"),
                   #legend.spacing.y = unit(0.75, "cm"),
                   #facets
                   strip.background = element_blank(),
                   strip.text = element_text(size = 8 ,face="bold", color = "black", hjust =0),
)



# Subset the data by species
macro <- subset(plot_dat, species == "macrocystis_pyrifera") %>%
              #fix names
        mutate(site = str_to_title(gsub("_", " ", site)),
         site = str_replace(site, "Dc", "DC"),
         site = str_replace(site, "Uc", "UC"))
strong <- subset(plot_dat, species == "strongylocentrotus_purpuratus")%>%
  #fix names
  mutate(site = str_to_title(gsub("_", " ", site)),
         site = str_replace(site, "Dc", "DC"),
         site = str_replace(site, "Uc", "UC"))


# Create the dumbbell plot with arrows
A <- ggplot(data = macro, aes(x = mean_count_after, y = mean_count_before)) +
  geom_segment(aes(x = strong$mean_count_after, xend = strong$mean_count_after, y = macro$mean_count_before, yend = macro$mean_count_after, color = mean_sim), size = 1, arrow = arrow(length = unit(0.25, "cm"))) +
  geom_point(aes(x = strong$mean_count_after, y = macro$mean_count_before, color = mean_sim), size = 3) +
  xlab("Mean purple sea urchin density 2014-2020 (per 60 m²)") +
  ylab("Kelp stipe density (per 60 m²)") +
  labs(color = "Community similarity \n (2007-2013 vs. 2014-2020)")+
  #scale_x_log10("Mean Count After (Strongylocentrotus Purpuratus)") +
  scale_color_gradient(name = "Mean Sim") +
  scale_color_viridis_c() +
  theme_classic()+
  #ggrepel::geom_label_repel(data = macro, aes(x = strong$mean_count_after, y = macro$mean_count_before, label = site),size=2, box.padding = 1, force = 50,min.segment.length = 1)+
  my_theme
A


##build schematic
# build schematic
B <- ggplot() +
  geom_segment(data = data.frame(x = c(0), y = c(2.7), xend = c(25), yend = c(2), facet_var = "i. no community structure change"),
               aes(x = x, y = y, xend = xend, yend = yend), color = "#008080", arrow = arrow(length = unit(0.25, "cm"), ends = "last"), size = 1) +
  geom_segment(data = data.frame(x = c(0), y = c(2.7), xend = c(25), yend = c(0.2), facet_var = "ii. kelp loss"),
               aes(x = x, y = y, xend = xend, yend = yend), color = "#FF7F00", arrow = arrow(length = unit(0.25, "cm"), ends = "last"), size = 1) +
  geom_point(data = data.frame(x = c(0), y = c(2.7), facet_var = "i. no community structure change"),
             aes(x = x, y = y), color = "#008080", size = 3) +
  geom_point(data = data.frame(x = c(0), y = c(2.7), facet_var = "ii. kelp loss"),
             aes(x = x, y = y), color = "#FF7F00", size = 3) +
  facet_wrap(~ facet_var, ncol = 2, labeller = labeller(facet_var = c("i. no community structure change" = "i. Kelp persistence \nsmall community structure change", "ii. kelp loss" = "ii. Kelp loss \nlarge community structure change"))) +
  xlab("Purple sea urchin density (per m²)") +
  ylab("Kelp stipe density (per m²)") +
  labs(tag = "A",
       title = "Hypothesized relationships")+
  labs(color = "Community similarity \n (2007-2013 vs. 2014-2020)") +
  scale_color_identity() +
  scale_x_continuous(limits = c(0, 30), breaks = seq(0, 30, 10)) +
  scale_y_continuous(limits = c(0, 3.5), breaks = seq(0, 3.5, 1)) +
  theme_classic() +
  my_theme + theme(aspect.ratio = 1.2,
                   plot.margin = margin(0, 50, 0, 0, "pt")) # decrease margin to align with C

B




# Create the dumbbell plot with arrows and log scale for strongylocentrotus_purpuratus
C <- ggplot(data = macro, aes(x = mean_count_after/60, y = mean_count_before/60)) +
  geom_segment(aes(x = strong$mean_count_before/60, xend = strong$mean_count_after/60, y = macro$mean_count_before/60, yend = macro$mean_count_after/60, color = mean_sim), size = 1, arrow = arrow(length = unit(0.25, "cm"))) +
  geom_point(aes(x = strong$mean_count_before/60, y = macro$mean_count_before/60, color = mean_sim), size = 3) +
  xlab("Purple sea urchin density (per m²)") +
  ylab("Kelp stipe density (per m²)") +
  labs(color = "Community similarity \n (2007-2013 vs. 2014-2020)",
       title = "Observed relationships") +
  scale_color_gradient2(low = "#FF7F00", high = "#008080", mid = "gray80", midpoint = 0.41) +  # Adjust the midpoint value
  theme_classic() +
  labs(tag = "B")+
  #ggrepel::geom_label_repel(data = macro, aes(x = strong$mean_count_before/60, y = macro$mean_count_before/60, label = site), size = 2, box.padding = 1, force = 20, min.segment.length = 1) +
  my_theme
C

C_no_legend <- C + theme(legend.position = "none")

# Get the legend from plot C
legend_C <- cowplot::get_legend(C)


D <- cowplot::plot_grid(B, C_no_legend, ncol = 1, nrow = 2, rel_heights = c(0.45,0.55), align = "hv",axis = 1) 

E <- cowplot::plot_grid(D, legend_C, ncol = 2, nrow = 1, align = "hv",axis = 1, rel_widths = c(0.5,0.1)) 
E



# Save the combined plot
ggsave(E, filename = file.path(figdir, "Fig3_dumbbell_new5.png"), 
       width = 7, height = 8, bg = "white", units = "in", dpi = 600)




