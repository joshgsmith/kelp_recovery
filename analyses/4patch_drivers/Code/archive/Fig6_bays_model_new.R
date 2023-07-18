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
#process response variable

response_vars <- swath_raw %>% dplyr::select(year, site, zone, transect,
                                             macrocystis_pyrifera) %>%
                  group_by(year, site) %>%
                  summarize(stipe_mean = mean(macrocystis_pyrifera, na.rm = TRUE))


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
#check for collinearity and normalize where needed

# Select the predictor variables to check for collinearity
#predictors <- c("beuti_avg", "cuti_avg", "sst_c_mean")

# Calculate the correlation matrix between the predictor variables
#cor_matrix <- cor(mod_dat_build3[predictors])

# Print the correlation matrix
#print(cor_matrix)

#normalize predictors to account for collinearity
#mod_dat_build3$sst_norm <- scale(mod_dat_build3$sst_c_mean)

################################################################################
#build model

# Specify the formula
formula <-  bf(stipe_mean ~ #resistance + 
                 (1 | site) + year + vrm_sum + bat_mean + beuti_month_obs +
                           npp_ann_mean + wave_hs_max + orb_vmax +
                           slope_mean + sst_month_obs + baseline_kelp + baseline_kelp_cv +
                           urchin_density)


# Standardize the predictors using z-scores
mod_dat_std <- mod_dat
mod_dat_std[, c("vrm_sum", "bat_mean", "beuti_month_obs", "npp_ann_mean",
                "wave_hs_max", "orb_vmax", "slope_mean", "sst_month_obs", 
                "baseline_kelp","urchin_density", "baseline_kelp_cv")] <- 
  scale(mod_dat[, c("vrm_sum", "bat_mean", "beuti_month_obs", "npp_ann_mean",
       "wave_hs_max", "orb_vmax", "slope_mean", "sst_month_obs", "baseline_kelp",
       "urchin_density","baseline_kelp_cv")])


################################################################################

# Define priors for regression coefficients
prior_coeff <- prior(normal(0, 5), class = b)

# Define prior for the standard deviation
prior_sd <- prior(cauchy(0, 2.5), class = sigma)


# Combine priors
priors <- c(prior_coeff, prior_sd)

set.seed(1985)
model <- brm(formula = formula,
             data = mod_dat_std,
             family = skew_normal(alpha = ~1),
             #family = gaussian(link = "log"), # to ensure non negative predictions
             prior = priors,
             warmup = 1000,
             iter = 2000,
             chains = 4,
             cores = 4)

# Run the model
fit <- update(model, iter = 4000)
samples <- posterior::as_draws(fit)

# Get the population-level effects for table


table_output <- capture.output(print(fit))


#writeLines(table_output, file.path(tab_dir, "TableS2_fit_output_table.txt"))



#posterior check
# Plot the posterior distributions for each parameter
#bayesplot::mcmc_areas(samples)

# Plot the posterior distributions for each parameter excluding "Intercept"
bayesplot::mcmc_areas(fit, pars = c("b_npp_ann_mean",
                                    "b_baseline_kelp",
                                    "b_baseline_kelp_cv",
                                    #"b_year",
                                    "b_vrm_sum",
                                    "b_bat_mean",
                                    "b_beuti_month_obs",
                                    #"b_npp_ann_mean",
                                    "b_wave_hs_max",
                                    "b_orb_vmax",
                                    "b_slope_mean",
                                    "b_sst_month_obs")) 

bayesplot::pp_check(fit, ndraws = 1000)

library(bayesplot)


y <- mod_dat_std$stipe_mean
yrep <- posterior_predict(fit, newdata = mod_dat_std, draws = 500, na.rm = TRUE)

ppc_stat(y, yrep, stat = "median")

################################################################################
#test rescale with beta 

library(brms)
library(scales)

mod_dat_std <- mod_dat
mod_dat_std[, c("vrm_sum", "bat_mean", "beuti_month_obs", "npp_ann_mean",
                "wave_hs_max", "orb_vmax", "slope_mean", "sst_month_obs", 
                "baseline_kelp","urchin_density", "baseline_kelp_cv")] <- 
  scale(mod_dat[, c("vrm_sum", "bat_mean", "beuti_month_obs", "npp_ann_mean",
                    "wave_hs_max", "orb_vmax", "slope_mean", "sst_month_obs", "baseline_kelp",
                    "urchin_density","baseline_kelp_cv")])

mod_dat_std$stipe_mean <- scales::rescale(mod_dat$stipe_mean, to = c(0.001, 0.99999))

# Specify the formula with rescaled stipe_mean
formula <- bf(stipe_mean ~ #resistance + 
                (1 | site) + year + 
                #(year | site) +
                vrm_sum + bat_mean + beuti_month_obs +
                npp_ann_mean + wave_hs_max + orb_vmax +
                slope_mean + sst_month_obs + baseline_kelp + baseline_kelp_cv +
                urchin_density)

# Calculate summary statistics from the data

b_mean <- mean(mod_dat_std$vrm_sum)
b_sd <- sd(mod_dat_std$vrm_sum)


# Define priors for regression coefficients
#prior_coeff <- prior(student_t(40, 0.5, 2.5), class = "Intercept")
prior_coeff <- prior(student_t(10, 0.24, 0.21), class = "Intercept")
prior_b <- prior(normal(0, 5), class = "b")
prior_phi <- prior(gamma(0.1, 0.1), class = "phi") #toy with this. Default is 0.1, 0.1

# Combine priors
priors <- c(prior_coeff, prior_b, prior_phi)

# Set the seed
set.seed(1985)

# Run the model
model <- brm(formula = formula,
             data = mod_dat_std,
             family = "beta",
             prior = priors,
             warmup = 1000,
             iter = 2000,
             chains = 4,
             cores = 4
             #control = list(adapt_delta = 0.8)
             )

# Update the model
fit <- update(model, iter = 8000)
samples <- posterior::as_draws(fit)

# Get the population-level effects for table
#summary(fit)

bayesplot::pp_check(fit, ndraws = 1000)
bayesplot::pp_check(fit, type = 'stat',stat='mean')

summary(fit)


################################################################################
########BEST FIT #1

###THE BEST THE BEST



library(brms)
library(scales)

mod_dat_std <- mod_dat
mod_dat_std[, c("vrm_sum", "bat_mean", "beuti_month_obs", "npp_ann_mean",
                "wave_hs_max", "orb_vmax", "slope_mean", "sst_month_obs", 
                "baseline_kelp","urchin_density", "baseline_kelp_cv")] <- 
  scale(mod_dat[, c("vrm_sum", "bat_mean", "beuti_month_obs", "npp_ann_mean",
                    "wave_hs_max", "orb_vmax", "slope_mean", "sst_month_obs", "baseline_kelp",
                    "urchin_density","baseline_kelp_cv")])

mod_dat_std <- mod_dat_std %>% filter(stipe_mean > 0)

#mod_dat_std$stipe_mean <- scales::rescale(mod_dat$stipe_mean, to = c(0, 0.9999))
hist(sqrt(mod_dat_std$stipe_mean))


# Specify the formula with rescaled stipe_mean
formula <- bf(sqrt(stipe_mean)  ~ #resistance + 
                (1 | site) + (1|year) + 
                #(year | site) +
                vrm_sum + bat_mean + beuti_month_obs +
                npp_ann_mean + wave_hs_max + orb_vmax +
                slope_mean + sst_month_obs + baseline_kelp + baseline_kelp_cv +
                urchin_density)


i_mean <- mean(sqrt(mod_dat_std$stipe_mean))
i_sd <- sd(sqrt(mod_dat_std$stipe_mean))

# Define priors for regression coefficients
#prior_coeff <- prior(student_t(40, 0.5, 2.5), class = "Intercept")
prior_coeff <- prior(student_t(5, 9, 4), class = "Intercept")
prior_b <- prior(normal(0, 5), class = "b")
#prior_phi <- prior(gamma(0.1, 0.1), class = "phi") #toy with this. Default is 0.1, 0.1

# Combine priors
priors <- c(prior_coeff, prior_b)

# Set the seed
set.seed(1985)

# Run the model
model <- brm(formula = formula,
             data = mod_dat_std,
             prior = priors,
            family = skew_normal(),
            #family = Beta(link = "logit", link_phi = "log")
           # negbinomial(link = "log", link_shape = "log"),
             warmup = 1000,
             iter = 10000,
             chains = 4,
             cores = 4
             #control = list(adapt_delta = 0.8)
)

# Update the model
fit <- update(model, iter = 5000)
samples <- posterior::as_draws(fit)

# Get the population-level effects for table
#summary(fit)

# Theme
my_theme <-  theme(axis.text=element_text(size=6, color = "black"),
                  # axis.text.y = element_blank(),
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

bayesplot::color_scheme_set(scheme = "blue") 
A <- bayesplot::pp_check(fit, ndraws = 1000) + 
  theme_bw() + 
  labs(tag = "A",
       title = "Posterior distribution, observed vs out-of-sample predictions")+
  my_theme  
A

B <- bayesplot::pp_check(fit, type = 'stat',stat='mean') + 
  labs(tag = "B",
       title = "Mean observed density vs. distribution of out-of-sample predictions")+
  theme_bw() + my_theme
B

################################################################################
#Plot

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
predictor_names <- c("b_npp_ann_mean" = "Net primary productivity", 
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
  scheme11 = c("#FF4081", "#F8BBD0", "#E91E63", "#FF80AB", "#FF80AB", "#EC407A")
)


# Function to plot posterior distribution for each parameter
plot_posterior <- function(parameter, color_scheme) {
  bayesplot::color_scheme_set(color_scheme)
  plot <- bayesplot::mcmc_areas(fit, pars = parameter) +
   #coord_cartesian(xlim = c(-1.5, 1)) +
    coord_cartesian(xlim = c(-4, 4)) +
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

kelp <- ggplot(data = mod_dat, aes(x = resistance, y = baseline_kelp/60)) +
  geom_boxplot(fill = "#D8B166", color = "black") +
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

sst <- ggplot(data = mod_dat, aes(x = resistance, y = sst_month_obs)) +
  geom_boxplot(fill = "#BA68C8", color = "black") +
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

beuti <- ggplot(data = mod_dat, aes(x = resistance, y = beuti_month_obs)) +
  geom_boxplot(fill = "#A5D6A7", color = "black") +
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

slope <- ggplot(data = mod_dat, aes(x = resistance, y = slope_mean)) +
  geom_boxplot(fill = "#E57373", color = "black") +
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

rugosity <- ggplot(data = mod_dat, aes(x = resistance, y = vrm_mean)) +
  geom_boxplot(fill = "#0D47A1", color = "black") +
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

orb_v <- ggplot(data = mod_dat, aes(x = resistance, y = orb_vmax)) +
  geom_boxplot(fill = "#E0E0E0", color = "black") +
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

bat <- ggplot(data = mod_dat, aes(x = resistance, y = bat_mean)) +
  geom_boxplot(fill = "#FFCC80", color = "black") +
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
  geom_boxplot(fill = "#81D4FA", color = "black") +
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


kelp_cv <- ggplot(data = mod_dat, aes(x = resistance, y = baseline_kelp_cv)) +
  geom_boxplot(fill = "#78909C", color = "black") +
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

wave_h <- ggplot(data = mod_dat, aes(x = resistance, y = wave_hs_max)) +
  geom_boxplot(fill = "#F48FB1", color = "black") +
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

urchin <- ggplot(data = mod_dat, aes(x = resistance, y = urchin_density/60)) +
  geom_boxplot(fill = "#FF4081", color = "black") +
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


predictors1 <- ggpubr::ggarrange(kelp, sst, beuti,slope,rugosity,  orb_v, bat, npp, 
                                  kelp_cv, wave_h, urchin, ncol=2, nrow=6, align = "v")  + 
  labs(tag = "B") + theme(plot.tag = element_text(size=8, face="plain")) 
predictors <- annotate_figure(predictors1,
                                           bottom = text_grob("Site type", 
                                                              hjust = 4.8, vjust = 0.1, x = 1, size = 10))
#predictors

full_plot <- ggarrange(g, predictors, nrow=1,  #widths=c(1.3,2)
                       widths = c(0.5,0.5),
                       heights = c(0.1,0.9)) + theme(plot.margin = margin(10, 1, 1, 1))
                       #plot.margin = margin(10, 10, 10, 10, "pt")
                       
full_plot



ggsave(full_plot, filename=file.path(figdir, "Fig5_predictors_new6.png"), 
       width=7, height=9, bg="white", units="in", dpi=600,
       device = "png")


