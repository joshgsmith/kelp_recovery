#Mixed effects model of patch drivers
#Joshua G. Smith
#May 9, 2023

rm(list=ls())
librarian::shelf(tidyverse, here, vegan, cowplot, ggpubr, lme4, rstan, brms)


################################################################################
#set directories and load data
basedir <- "/Volumes/seaotterdb$/kelp_recovery/"
figdir <- here::here("analyses","4patch_drivers","Figures")
outdir <- here::here("analyses","4patch_drivers","Output")
tab_dir <- here::here("analyses","4patch_drivers","Tables")

#load sampling data
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

#load model predictors
mod_predict <- readRDS(file.path(basedir, "data/environmental_data/predictors_at_pisco_sites_v2.Rds"))

#load published urchin behavior data
#accessed from "https://github.com/joshgsmith/PatchDynamics/blob/main/Data/raw_data/patch_data.csv"
patch_raw <- read.csv(file.path(outdir, "patch_data.csv"))


################################################################################
#process environmental data

mod_predict_build1 <- mod_predict %>%
  dplyr::filter(year>=2007 & year <= 2020)%>%
  #calculate annual means
  group_by(site, year)%>%
  summarize(across(8:46,mean, na.rm=TRUE)) %>%
  #drop monthly statistics
  dplyr::select(!(c(cuti_month_baseline, cuti_month_sd,
                    beuti_month_baseline, beuti_month_baseline_sd,
                    sst_month_baseline,
                    sst_month_baseline_sd, sst_month_anom, 
  ))) %>% #these are calculated at monthly intervals so irrelevant here. 
  dplyr::mutate_all(~ ifelse(is.nan(.), NA, .))

################################################################################
#calculate baseline kelp density

kelp_baseline <- swath_raw %>% dplyr::select(year, site, zone, transect, macrocystis_pyrifera)%>%
  filter(year <= 2013)%>%
  group_by(site)%>%
  summarize(baseline_kelp = mean(macrocystis_pyrifera, na.rm=TRUE),
            baseline_kelp_cv = (sd(macrocystis_pyrifera, na.rm = TRUE) / mean(macrocystis_pyrifera, na.rm = TRUE)) * 100
  )

################################################################################
#calculate mean urchin densities

urchin_density <- swath_raw %>% dplyr::select(year, site, zone, transect, strongylocentrotus_purpuratus)%>%
  group_by(year, site)%>%
  summarize(urchin_density = mean(strongylocentrotus_purpuratus, na.rm=TRUE))

################################################################################
#calculated inferred behavior

##Note: prop_pur_exp is average across replicate transects, not calculated directly 
#from the density averages. 

#examine data

ggplot(data = patch_raw, aes(x = mac_stipe_den, y = prop_pur_exp)) +
  geom_point() +
  xlab("mac_stipe_den") +
  ylab("prop_pur_exposed") +
  ggtitle("Scatter Plot of prop_pur_exposed vs. mac_stipe_den")

#remove six outliers

exp_build1 <- patch_raw %>% filter(!(prop_pur_exp > 0.25 & mac_stipe_den > 4))

# Remove rows with NA values in mac_stipe_den column
exp_build1 <- exp_build1[!is.na(exp_build1$mac_stipe_den), ]

#now fit model
# Fit negative exponential curve using nls()
#fit <- nls(prop_pur_exp ~ a * exp(-b * mac_stipe_den), data = exp_build1, start = list(a = 1, b = 1))

# Fit modified exponential decay model using nls()
fit <- nls(prop_pur_exp ~ a * (1 - exp(-b * mac_stipe_den)) + c, data = exp_build1, start = list(a = 1, b = 1, c = 0))

# Create line data using the cleaned exp_build1 data frame
line_data <- data.frame(mac_stipe_den = seq(min(exp_build1$mac_stipe_den), max(exp_build1$mac_stipe_den), length.out = 100))
line_data$prop_pur_exposed <- predict(fit, newdata = line_data)

# Plot with fitted line
ggplot(data = exp_build1, aes(x = mac_stipe_den, y = prop_pur_exp)) +
  geom_point() +
  geom_line(data = line_data, aes(x = mac_stipe_den, y = prop_pur_exposed), color = "red") +
  xlab("mac_stipe_den") +
  ylab("prop_pur_exposed") +
  ggtitle("Scatter Plot of prop_pur_exposed vs. mac_stipe_den")

# Calculate R-squared value
residuals <- residuals(fit)
y_mean <- mean(exp_build1$prop_pur_exp, na.rm=TRUE)
ss_total <- sum((exp_build1$prop_pur_exp - y_mean)^2, na.rm=TRUE)
ss_residual <- sum(residuals^2)
r_squared <- 1 - (ss_residual / ss_total)

# Get the formula
formula_text <- as.character(as.formula(fit))


################################################################################
#process response variable and exposed urchin 

response_vars <- swath_raw %>% dplyr::select(year, site, zone, transect,
                                             macrocystis_pyrifera) %>%
  group_by(year, site) %>%
  summarize(stipe_mean = mean(macrocystis_pyrifera, na.rm = TRUE),
            prop_urch_exp = -0.6487 * (1 - exp(-5.1438 * (stipe_mean / 60))) + 0.7169)


################################################################################
#join sampling data with environmental

mod_dat <- left_join(response_vars, mod_predict_build1, by=c("year","site")) %>%
  #create lagged variables
  mutate(#npp_lag1 = lag(npp,1),
    #npp_lag2 = lag(npp,2),
    #designate resistant site
    resistance = ifelse(site == "HOPKINS_UC" |
                          site == "CANNERY_UC" |
                          site == "SIREN" |
                          site == "CANNERY_DC",
                        #site == "BUTTERFLY_DC",
                        "resistant","transitioned")
  ) %>%
  left_join(kelp_baseline, by="site") %>% 
  left_join(urchin_density, by = c("year","site"))

################################################################################
#prep space

# Use the cmdstanr backend for Stan 
options(mc.cores = 6,
        brms.backend = "cmdstanr")

# Set some global Stan options
CHAINS <- 4
ITER <- 10000
WARMUP <- 2000
BAYES_SEED <- 1985

# Use the Johnson color palette
clrs <- MetBrewer::met.brewer("Hokusai3")

# Tell bayesplot to use the Johnson palette (for things like pp_check())
bayesplot::color_scheme_set(c(clrs[1], clrs[2], clrs[3], clrs[4], clrs[5],clrs[6]))


################################################################################
#Build model


mod_dat_std <- mod_dat 
mod_dat_std[, c("prop_urch_exp","vrm_sum", "bat_mean", "beuti_month_obs", "npp_ann_mean",
                "wave_hs_max", "orb_vmax", "slope_mean", "sst_month_obs", 
                "baseline_kelp","urchin_density", "baseline_kelp_cv")] <- 
  scale(mod_dat[, c("prop_urch_exp","vrm_sum", "bat_mean", "beuti_month_obs", "npp_ann_mean",
                    "wave_hs_max", "orb_vmax", "slope_mean", "sst_month_obs", "baseline_kelp",
                    "urchin_density","baseline_kelp_cv")])

mod_dat_std$stipe_mean <-round(mod_dat_std$stipe_mean,0)

hist(mod_dat_std$stipe_mean)

hurdle_poisson_full <- brm(
  bf(stipe_mean  ~ 
       vrm_sum + #bat_mean + 
       beuti_month_obs +
       #npp_ann_mean + 
       wave_hs_max + orb_vmax +
        slope_mean + sst_month_obs + baseline_kelp + baseline_kelp_cv +
       urchin_density + 
       #year + (1 | year / site), 
        year + (1|year / site),
     hu ~  vrm_sum + #bat_mean + 
       beuti_month_obs +
       #npp_ann_mean + 
       wave_hs_max + orb_vmax +
       slope_mean + 
       sst_month_obs + 
       baseline_kelp + baseline_kelp_cv +
       urchin_density + 
       #year + (1 | year / site)
       year + (1|year / site)
  ),
  data = mod_dat_std,
  family = hurdle_poisson(),
  #control = list(adapt_delta = 0.9),
  chains = CHAINS, iter = ITER, warmup = WARMUP, seed = BAYES_SEED
)

hurdle_poisson_red1 <- brm(
  bf(stipe_mean  ~ 
       vrm_sum + #bat_mean + 
       beuti_month_obs +
       #npp_ann_mean + 
       wave_hs_max + #orb_vmax +
       slope_mean + sst_month_obs + baseline_kelp + #baseline_kelp_cv +
       urchin_density + 
       #year + (1 | year / site), 
       year + (1|year / site),
     hu ~  vrm_sum + #bat_mean + 
       beuti_month_obs +
       #npp_ann_mean + 
       wave_hs_max + #orb_vmax +
       slope_mean + 
       sst_month_obs + 
       baseline_kelp + #baseline_kelp_cv +
       urchin_density + 
       #year + (1 | year / site)
       year + (1|year / site)
  ),
  data = mod_dat_std,
  family = hurdle_poisson(),
  #control = list(adapt_delta = 0.9),
  chains = CHAINS, iter = ITER, warmup = WARMUP, seed = BAYES_SEED
)

hurdle_poisson_red2 <- brm(
  bf(stipe_mean  ~ 
       prop_urch_exp +
       vrm_sum + bat_mean + 
       beuti_month_obs +
       npp_ann_mean + 
       wave_hs_max + orb_vmax +
       slope_mean + sst_month_obs + baseline_kelp + baseline_kelp_cv +
       urchin_density + 
       year + (1|year / site),
     hu ~  #vrm_sum + #bat_mean + 
       prop_urch_exp +
       #beuti_month_obs +
       #npp_ann_mean + 
       #wave_hs_max + #orb_vmax +
       #slope_mean + 
       #sst_month_obs + 
       #baseline_kelp + #baseline_kelp_cv +
       urchin_density + 
       year + (1|year / site)
  ),
  data = mod_dat_std,
  family = hurdle_poisson(),
  control = list(adapt_delta = 0.9),
  chains = CHAINS, iter = ITER, warmup = WARMUP, seed = BAYES_SEED
)

# posterior checks
# Exponential
pp_check(hurdle_poisson_red2, ndraws = 100)

pred <- posterior_predict(hurdle_poisson_red2)
bayesplot::ppc_dens_overlay(y = log1p(mod_dat_std$stipe_mean), 
                            yrep = log1p(pred[1:100,]))

bayesplot::ppc_dens_overlay(y = mod_dat_std$stipe_mean, 
                            yrep = pred[1:100,])


my_theme <-  theme(axis.text=element_text(size=6, color = "black"),
                   #axis.text.y = element_blank(),
                   axis.title=element_text(size=8,color = "black"),
                   plot.tag=element_text(size=8, face = "plain",color = "black"),
                   plot.title =element_text(size=7, face="plain",color = "black", vjust=-1),
                   # Gridlines 
                   panel.grid.major = element_blank(), 
                   panel.grid.minor = element_blank(),
                   panel.background = element_blank(), 
                   axis.line = element_line(colour = "black"),
                   # Legend
                   legend.key = element_blank(),
                   legend.background = element_rect(fill=alpha('blue', 0)),
                   legend.key.height = unit(1, "lines"), 
                   legend.text = element_text(size = 6,color = "black"),
                   legend.title = element_text(size = 7,color = "black"),
                   #legend.spacing.y = unit(0.75, "cm"),
                   #facets
                   strip.background = element_blank(),
                   strip.text = element_text(size = 6 ,face="plain",color = "black"),
)

bayesplot::color_scheme_set(scheme = "blue") 
A <- bayesplot::pp_check(hurdle_poisson_red2, ndraws = 1000) + 
  theme_bw() + 
  labs(tag = "A",
       title = "Posterior distribution, observed vs out-of-sample predictions")+
  xlab("Density of kelp stipes (per 60 m² transect)")+
  ylab("Frequency")+
  scale_x_continuous(limits = c(0, 400))+
  my_theme  

B <- bayesplot::pp_check(hurdle_poisson_red2, type = 'stat',stat='mean') + 
  labs(tag = "B",
       title = "Mean observed density vs. distribution of out-of-sample predictions")+
  xlab("Density of kelp stipes (per 60 m² transect)")+
  ylab("Frequency")+
  theme_bw() + my_theme

g <- ggpubr::ggarrange(A,B,ncol = 1)
g

ggsave(g, filename=file.path(figdir, "FigS4_pp_checks.png"), 
       width=7, height=9, bg="white", units="in", dpi=600,
       device = "png")

table_output <- capture.output(print(hurdle_poisson_red2))
writeLines(table_output, file.path(tab_dir, "TableS2_fit_output_table.txt"))



################################################################################
#plot

fit <- hurdle_poisson_red2 

# Theme
my_theme <-  theme(axis.text=element_text(size=6, color = "black"),
                   axis.text.y = element_blank(),
                   axis.title=element_text(size=8,color = "black"),
                   plot.tag=element_text(size=8, face = "bold",color = "black"),
                   plot.title =element_text(size=7, face="bold",color = "black", vjust=-1),
                   # Gridlines 
                   panel.grid.major = element_blank(), 
                   panel.grid.minor = element_blank(),
                   panel.background = element_blank(), 
                   axis.line = element_line(colour = "black"),
                   # Legend
                   legend.key = element_blank(),
                   legend.background = element_rect(fill=alpha('blue', 0)),
                   legend.key.height = unit(1, "lines"), 
                   legend.text = element_text(size = 6,color = "black"),
                   legend.title = element_text(size = 7,color = "black"),
                   #legend.spacing.y = unit(0.75, "cm"),
                   #facets
                   strip.background = element_blank(),
                   strip.text = element_text(size = 6 ,face="plain",color = "black"),
)



# Map predictor names
predictor_names <- c("b_prop_urch_exp" = "Behavior (proportion active urchins)", 
                     "b_npp_ann_mean" = "Net primary productivity", 
                     "b_baseline_kelp" = "Baseline kelp density", 
                     "b_baseline_kelp_cv" = "Baseline kelp coefficient of variation",
                     "b_urchin_density" = "Urchin density",
                     #"b_year",
                     "b_vrm_sum" = "Rugosity", 
                     "b_bat_mean" = "Mean depth (m)", 
                     "b_beuti_month_obs" = "Upwelling (BEUTI)",
                     # "b_npp_ann_mean" = "", 
                     "b_wave_hs_max" = "Wave height (m)", 
                     "b_orb_vmax" = "Orbital velocity",
                     "b_slope_mean" = "Reef slope", 
                     "b_sst_month_obs" = "Sea surface temperature (°C)")

# Extract the posterior samples
posterior_samples <- as.matrix(fit, pars = names(predictor_names))

# Calculate the mean for each parameter
means <- colMeans(posterior_samples)

# Sort the parameter names based on mean values in descending order
sorted_pars <- names(means)[order(-means)]

# Define the number of color variations for each parameter
num_variations <- 6 #9C27B0

# Define the manually created color schemes
color_schemes <- list(
  scheme1 = c("#D8B166", "#FFECB3", "#E6A800", "#FFD633", "#FFD633", "#E6C300"),
  scheme2 = c("#BA68C8", "#E1BEE7", "#9C27B0", "#8E24AA", "#8E24AA", "#CE93D8"),
  scheme3 = c("#A5D6A7", "#C8E6C9", "#4CAF50", "#2E7D32", "#2E7D32", "#81C784"),
  scheme4 = c("#E57373", "#FFCDD2", "#F44336", "#C62828", "#C62828", "#EF9A9A"),
  scheme5 = c("#0D47A1", "#0D47A1", "#1976D2", "#1976D2", "#1976D2", "#64B5F6"),
  scheme6 = c("#EEEEEE", "#F5F5F5", "#9E9E9E", "#757575", "#757575", "#E0E0E0"),
  scheme7 = c("#FFCC80", "#FFE0B2", "#FF9800", "#E65100", "#E65100", "#FFB74D"),
  scheme8 = c("#81D4FA", "#B3E5FC", "#03A9F4", "#0288D1", "#0288D1", "#4FC3F7"),
  scheme9 = c("#78909C", "#CFD8DC", "#607D8B", "#455A64", "#455A64", "#78909C"),
  scheme10 = c("#F48FB1", "#F8BBD0", "#E91E63", "#C2185B", "#C2185B", "#F06292"),
  scheme11 = c("#FF4081", "#F8BBD0", "#E91E63", "#FF80AB", "#FF80AB", "#EC407A"),
  scheme12 = c("#5C6BC0", "#D1C4E9", "#8E24AA", "#9C27B0", "#BA68C8", "#BA68C8")
)


# Function to plot posterior distribution for each parameter
plot_posterior <- function(parameter, color_scheme) {
  bayesplot::color_scheme_set(color_scheme)
  plot <- bayesplot::mcmc_areas(fit, pars = parameter) +
    #coord_cartesian(xlim = c(-1.5, 1)) +
    coord_cartesian(xlim = c(-2.2, 1)) +
    theme(plot.margin = margin(1, 10, 5, 10)) +
    labs(y = NULL) +
    labs(title = predictor_names[parameter]) +
    my_theme +
    theme(axis.text.y = element_blank()) +
    theme_bw() +
    my_theme
  
  # Add a darker vertical line at 0
  plot <- plot + geom_vline(xintercept = 0, color = "navyblue", linetype = "solid", size = 1) +
    theme(plot.margin = margin(1, 10, 1, 10)) 
  
  return(plot)
}

# Generate posterior plots for each parameter
plot_list <- list()
for (i in seq_along(sorted_pars)) {
  color_scheme <- color_schemes[[paste0("scheme", i)]]
  plot_list[[i]] <- plot_posterior(sorted_pars[i], color_scheme)
}

# Arrange the plots in a grid
g1 <- ggpubr::ggarrange(plotlist = plot_list, ncol = 1) + labs(tag = "A") + theme(plot.tag = element_text(size = 8, face = "plain"))
g <- ggpubr::annotate_figure(g1, left = text_grob("Density", 
                                                  rot = 90, vjust = 1, hjust=0.3, size = 10),
                             bottom = text_grob("Standardized beta coefficient", hjust=0.5, vjust=0, size = 10))
g


###############


# Theme 2
my_theme <-  theme(axis.text=element_text(size=6,color = "black"),
                   #axis.text.y = element_blank(),
                   axis.title=element_text(size=8,color = "black"),
                   plot.tag=element_text(size=8, face = "bold",color = "black"),
                   plot.title =element_text(size=7, face="bold",color = "black"),
                   # Gridlines 
                   panel.grid.major = element_blank(), 
                   panel.grid.minor = element_blank(),
                   panel.background = element_blank(), 
                   axis.line = element_line(colour = "black"),
                   # Legend
                   legend.key = element_blank(),
                   #legend.background = element_rect(fill=alpha('blue', 0)),
                   legend.key.height = unit(1, "lines"), 
                   legend.text = element_text(size = 6,color = "black"),
                   legend.title = element_text(size = 7,color = "black"),
                   #legend.spacing.y = unit(0.75, "cm"),
                   #facets
                   strip.background = element_blank(),
                   strip.text = element_text(size = 6 ,face="plain",color = "black"),
                   plot.margin = margin(0.2,1,0.2,1)
)

kelp <- ggplot(data = mod_dat %>% 
                 #take mean across since, since this is static
                 group_by(site, resistance) %>% summarise(baseline_kelp = mean(baseline_kelp)), 
               aes(x = resistance, y = baseline_kelp/60)) +
  geom_boxplot(fill = color_schemes$scheme1[1], color = "black") +
  geom_jitter(width = 0.1, height = 0.3, alpha = 0.2, size=1) +
  ggsignif::geom_signif(comparisons = list(c("resistant", "transitioned")),
                        map_signif_level = TRUE,
                        tip_length = c(0.01, 0.01),
                        textsize=3)+
  ylim(1,4)+
  xlab("") +
  ylab("Kelp baseline density \n(no. stipes per m²)") +
  ggtitle("Baseline kelp \ndensity") +
  labs(tag="")+
  theme_classic()+
  my_theme+
  scale_x_discrete(labels = c("Persistent", "Transitioned"))   # Renaming levels
#kelp

orb_v <- ggplot(data = mod_dat, aes(x = resistance, y = orb_vmax)) +
  geom_boxplot(fill = color_schemes$scheme2[1], color = "black") +
  geom_jitter(width = 0.1, height = 0.3, alpha = 0.2, size=1) +
  ggsignif::geom_signif(comparisons = list(c("resistant", "transitioned")),
                        map_signif_level = TRUE,
                        tip_length = c(0.01, 0.01),
                        textsize=3)+
  #ylim(0,20)+
  xlab("") +
  ylab("Orbital \nvelocity") +
  ggtitle("Orbital \nvelocity") +
  ylim(0,2500)+
  #labs(tag="B")+
  theme_classic()+
  my_theme+
  scale_x_discrete(labels = c("Persistent", "Transitioned"))  # Renaming levels
#orb_v

sst <- ggplot(data = mod_dat, aes(x = resistance, y = sst_month_obs)) +
  geom_boxplot(fill = color_schemes$scheme3[1], color = "black") +
  geom_jitter(width = 0.1, height = 0.3, alpha = 0.2, size=1) +
  ggsignif::geom_signif(comparisons = list(c("resistant", "transitioned")),
                        map_signif_level = TRUE,
                        y_position = 15.3,
                        tip_length = c(0.01, 0.01),
                        textsize=3)+
  ylim(12,16)+
  xlab("") +
  ylab("Sea surface temperature \n(°C, annual mean)") +
  ggtitle("Sea surface \ntemperature (°C)") +
  labs(tag="")+
  theme_classic()+
  my_theme+
  scale_x_discrete(labels = c("Persistent", "Transitioned")) # Renaming levels
#sst

slope <- ggplot(data = mod_dat %>%
                  #take mean across since, since this is static
                  group_by(site, resistance) %>% summarise(slope_mean = mean(slope_mean)), 
                aes(x = resistance, y = slope_mean)) +
  geom_boxplot(fill = color_schemes$scheme4[1], color = "black") +
  geom_jitter(width = 0.1, height = 0.3, alpha = 0.2, size=1) +
  ggsignif::geom_signif(comparisons = list(c("resistant", "transitioned")),
                        map_signif_level = TRUE,
                        tip_length = c(0.01, 0.01),
                        textsize=3)+
  ylim(0,20)+
  xlab("") +
  ylab("Slope \n(mean per site)") +
  ggtitle("Reef \nslope") +
  #labs(tag="B")+
  theme_classic()+
  my_theme+
  scale_x_discrete(labels = c("Persistent", "Transitioned"))  # Renaming levels
#slope

rugosity <- ggplot(data = mod_dat %>%
                     #take mean across since, since this is static
                     group_by(site, resistance) %>% summarise(vrm_mean = mean(vrm_mean)), 
                   aes(x = resistance, y = vrm_mean)) +
  geom_boxplot(fill = color_schemes$scheme5[1], color = "black") +
  geom_jitter(width = 0.1, height = 0, alpha = 0.2, size=1) +
  ggsignif::geom_signif(comparisons = list(c("resistant", "transitioned")),
                        map_signif_level = TRUE,
                        tip_length = c(0.01, 0.01),
                        textsize=3)+
  #ylim(50,220)+
  xlab("") +
  ylim(1.0,1.12)+
  ylab("Vector \nruggedness") +
  ggtitle("Rugosity") +
  labs(tag="")+
  theme_classic()+
  my_theme+
  scale_x_discrete(labels = c("Persistent", "Transitioned"))  # Renaming levels
#rugosity


beuti <- ggplot(data = mod_dat, aes(x = resistance, y = beuti_month_obs)) +
  geom_boxplot(fill = color_schemes$scheme6[1], color = "black") +
  geom_jitter(width = 0.1, height = 0.3, alpha = 0.2, size=1) +
  ggsignif::geom_signif(comparisons = list(c("resistant", "transitioned")),
                        map_signif_level = TRUE,
                        y_position = 11,
                        tip_length = c(0.01, 0.01),
                        textsize=3)+
  ylim(2,13)+
  xlab("") +
  ylab("BEUTI \n(annual mean") +
  ggtitle("Upwelling \n(BEUTI)") +
  labs(tag="")+
  theme_classic()+
  my_theme+
  scale_x_discrete(labels = c("Persistent", "Transitioned")) # Renaming levels
#beuti

kelp_cv <- ggplot(data = mod_dat %>%
                    #take mean across since, since this is static
                    group_by(site, resistance) %>% summarise(baseline_kelp_cv = mean(baseline_kelp_cv)), 
                  aes(x = resistance, y = baseline_kelp_cv)) +
  geom_boxplot(fill = color_schemes$scheme7[1], color = "black") +
  geom_jitter(width = 0.1, height = 0.3, alpha = 0.2, size=1) +
  ggsignif::geom_signif(comparisons = list(c("resistant", "transitioned")),
                        map_signif_level = TRUE,
                        tip_length = c(0.01, 0.01),
                        textsize=3)+
  #ylim(0,20)+
  xlab("") +
  ylim(40,175)+
  ylab("Coefficient of var.") +
  ggtitle("Baseline kelp \ncoefficient of var.") +
  #labs(tag="B")+
  theme_classic()+
  my_theme+
  scale_x_discrete(labels = c("Persistent", "Transitioned"))  # Renaming levels
#kelp_cv

bat <- ggplot(data = mod_dat %>%
                #take mean across since, since this is static
                group_by(site, resistance) %>% summarise(bat_mean = mean(bat_mean)), 
              aes(x = resistance, y = bat_mean)) +
  geom_boxplot(fill = color_schemes$scheme8[1], color = "black") +
  geom_jitter(width = 0.1, height = 0.3, alpha = 0.2, size=1) +
  ggsignif::geom_signif(comparisons = list(c("resistant", "transitioned")),
                        map_signif_level = TRUE,
                        tip_length = c(0.01, 0.01),
                        textsize=3)+
  ylim(5,23)+
  xlab("") +
  ylab("Depth (m)") +
  ggtitle("Mean \ndepth (m)") +
  labs(tag="")+
  theme_classic()+
  my_theme+
  scale_x_discrete(labels = c("Persistent", "Transitioned"))  # Renaming levels
#bat

npp <- ggplot(data = mod_dat, aes(x = resistance, y = npp_ann_mean)) +
  geom_boxplot(fill = color_schemes$scheme9[1], color = "black") +
  geom_jitter(width = 0.1, height = 0.3, alpha = 0.2, size=1) +
  ggsignif::geom_signif(comparisons = list(c("resistant", "transitioned")),
                        map_signif_level = TRUE,
                        y_position = 4000,
                        tip_length = c(0.01, 0.01),
                        textsize=3)+
  ylim(0,4500)+
  xlab("") +
  ylab("Net primary productivity \n(annual mean)") +
  ggtitle("Net primary \nproductivity") +
  labs(tag="")+
  theme_classic()+
  my_theme+
  scale_x_discrete(labels = c("Persistent", "Transitioned"))
#npp


urchin <- ggplot(data = mod_dat, aes(x = resistance, y = urchin_density/60)) +
  geom_boxplot(fill = color_schemes$scheme10[1], color = "black") +
  geom_jitter(width = 0.1, height = 0.3, alpha = 0.2, size=1) +
  ggsignif::geom_signif(comparisons = list(c("resistant", "transitioned")),
                        map_signif_level = TRUE,
                        tip_length = c(0.01, 0.01),
                        textsize=3)+
  ylim(0,70)+
  xlab("") +
  ylab("Urchin density \n(no. per m²)") +
  ggtitle("Urchin \ndensity") +
  #labs(tag="B")+
  theme_classic()+
  my_theme+
  scale_x_discrete(labels = c("Persistent", "Transitioned"))  # Renaming levels
#urchin



wave_h <- ggplot(data = mod_dat, aes(x = resistance, y = wave_hs_max)) +
  geom_boxplot(fill = color_schemes$scheme11[1], color = "black") +
  geom_jitter(width = 0.1, height = 0.3, alpha = 0.2, size=1) +
  ggsignif::geom_signif(comparisons = list(c("resistant", "transitioned")),
                        map_signif_level = TRUE,
                        tip_length = c(0.01, 0.01),
                        textsize=3)+
  #ylim(0,20)+
  xlab("") +
  ylim(0,13)+
  ylab("Wave \nheight (m)") +
  ggtitle("Wave \nheight") +
  #labs(tag="B")+
  theme_classic()+
  my_theme+
  scale_x_discrete(labels = c("Persistent", "Transitioned"))  # Renaming levels
#wave_h


behavior <- ggplot(data = mod_dat, aes(x = resistance, y = prop_urch_exp)) +
  geom_boxplot(fill = color_schemes$scheme12[1], color = "black") +
  geom_jitter(width = 0.1, height = 0.05, alpha = 0.2, size=1) +
  ggsignif::geom_signif(comparisons = list(c("resistant", "transitioned")),
                        map_signif_level = TRUE,
                        tip_length = c(0.01, 0.01),
                        textsize=3)+
  #ylim(0,20)+
  xlab("") +
  ylim(0,1.1)+
  ylab("Proportion active \nsea urchins") +
  ggtitle("Behavior") +
  #labs(tag="B")+
  theme_classic()+
  my_theme+
  scale_x_discrete(labels = c("Persistent", "Transitioned"))  # Renaming levels
#behavior



predictors1 <- ggpubr::ggarrange(kelp,orb_v, sst, slope, rugosity, beuti, kelp_cv, bat, npp, 
                                   urchin, wave_h, behavior, ncol=2, nrow=6, align = "v")  + 
  labs(tag = "B") + theme(plot.tag = element_text(size=8, face="plain")) 
predictors <- annotate_figure(predictors1,
                              bottom = text_grob("Site type", 
                                                 hjust = 3.2, vjust = 0.1, x = 1, size = 10))
#predictors

full_plot <- ggarrange(g, predictors, nrow=1,  #widths=c(1.3,2)
                       widths = c(0.5,0.5),
                       heights = c(0.1,0.9)) + theme(plot.margin = margin(10, 1, 1, 1))
#plot.margin = margin(10, 10, 10, 10, "pt")

full_plot



ggsave(full_plot, filename=file.path(figdir, "Fig5_predictors_new8.png"), 
       width=7, height=9, bg="white", units="in", dpi=600,
       device = "png")




