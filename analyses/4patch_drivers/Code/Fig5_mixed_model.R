#Mixed effects model of patch drivers
#Joshua G. Smith
#May 9, 2023

rm(list=ls())
librarian::shelf(tidyverse, here, vegan, cowplot, ggpubr, lme4)

library(dplyr)

################################################################################
#set directories and load data
basedir <- "/Volumes/seaotterdb$/kelp_recovery/"
figdir <- here::here("analyses","4patch_drivers","Figures")

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

mod_predict <- readRDS(file.path(basedir, "data/environmental_data/predictors_at_pisco_sites.Rds"))


################################################################################
#process environmental data

mod_predict_build1 <- mod_predict %>%
                        dplyr::filter(year>=2007 & year <= 2020)%>%
                        #calculate annual means
                        group_by(site, year)%>%
                        summarize(across(8:36,mean, na.rm=TRUE)) %>%
                        #drop monthly statistics
                        dplyr::select(!(c(cuti_month_baseline, cuti_month_sd,
                                          beuti_month_baseline, beuti_month_baseline_sd,
                                          sst_month_baseline,
                                          sst_month_baseline_sd, sst_month_anom, 
                                          ))) #these are calculated at monthly intervals so irrelevant here. 

################################################################################
#calculate baseline kelp density

kelp_baseline <- swath_raw %>% dplyr::select(year, site, zone, transect, macrocystis_pyrifera)%>%
                  filter(year <= 2013)%>%
                  group_by(site)%>%
                  summarize(baseline_kelp = mean(macrocystis_pyrifera, na.rm=TRUE))

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
                  mutate(npp_lag1 = lag(npp,1),
                         npp_lag2 = lag(npp,2),
                  #designate resistant site
                        resistance = ifelse(site == "HOPKINS_UC" |
                                              site == "CANNERY_UC" |
                                              site == "SIREN" |
                                              site == "CANNERY_DC",
                                              #site == "BUTTERFLY_DC",
                                            "resistant","transitioned")
                  ) %>%
          left_join(kelp_baseline, by="site")


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
#build full model

# Fit mixed model with random intercepts and slopes for site and year
full_mod <- lme4::lmer(stipe_mean ~ npp + vrm_sum + bat_mean + beuti_month_obs + 
                         slope_mean + sst_month_obs + year + baseline_kelp +
                         (1|site), data = mod_dat)

# fit the model
fit <- summary(full_mod)

# Print model summary
summary(full_mod)


################################################################################
#build sub models 

sub_mod1 <- lme4::lmer(stipe_mean ~year +
                         baseline_kelp +
                         beuti_month_obs + 
                         vrm_sum +
                         (1 | site), data = mod_dat)


# Compare the full model to the sub-models using an LRT
anova(full_mod, sub_mod1)




################################################################################
#prep plotting


# Extract the estimated effect size and confidence intervals for each predictor
effect_sizes <- data.frame(coef(summary(full_mod))[,1:3])
effect_sizes$predictor <- rownames(effect_sizes)

# Calculate confidence intervals
pred_ci <- confint(full_mod, level = 0.95) %>% data.frame() %>% rownames_to_column()


mod_out <- left_join(effect_sizes, pred_ci, by=c("predictor"="rowname")) %>%
              rename("ci_lower" = `X2.5..`,
                     "ci_upper" = `X97.5..`)


# Calculate the upper and lower bounds of the error bars
mod_out$lower_error <- effect_sizes$Estimate - 1.96 * effect_sizes$Std..Error
mod_out$upper_error <- effect_sizes$Estimate + 1.96 * effect_sizes$Std..Error

# Exclude the row corresponding to the intercept
mod_out <- mod_out[mod_out$predictor != "(Intercept)",]

################################################################################
#Plot

# Theme
my_theme <-  theme(axis.text=element_text(size=6),
                   axis.text.y = element_text(angle = 90, hjust = 0.5),
                   axis.title=element_text(size=8),
                   plot.tag=element_text(size=8),
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

# Create a forest plot with points and error bars
p1 <- ggplot(mod_out %>%
               
               mutate(predictor = ifelse(predictor == "sst_month_obs","SST",
                                         ifelse(predictor == "npp_lag2","NPP",
                                                ifelse(predictor == "beuti_month_obs","BEUTI",
                                                       ifelse(predictor == "baseline_kelp","Kelp baseline \n(2007-2013",
                                                              ifelse(predictor == "vrm_sum","Rugosity (m)",
                                                                     ifelse(predictor == "slope_mean","Slope (m)",
                                                                            ifelse(predictor == "bat_mean","Depth range (m)",
                                                                                   ifelse(predictor == "year","Year",predictor)))))))))
             
             , aes(x = reorder(predictor, -Estimate), y = Estimate)) +
  geom_point(#aes(color = Estimate), 
             size = 2) +
  #scale_color_gradient2(low = "navy", mid = "blue", high = "darkred", midpoint = 0) +
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0.2, size = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  coord_flip() +
  xlab("Predictor") +
  ylab("Effect size") +
  labs(color = "Effect size", tag = "A")+
  ggtitle("Mixed model effect sizes") +
  theme_classic() +
  my_theme

p1


slope <- ggplot(data = mod_dat, aes(x = resistance, y = slope_mean)) +
  geom_boxplot(fill = "#1B9E77", color = "black") +
  geom_jitter(width = 0.1, height = 0.3, alpha = 0.2, size=1) +
  ggsignif::geom_signif(comparisons = list(c("resistant", "transitioned")),
              map_signif_level = TRUE,
              tip_length = c(0.01, 0.01),
              textsize=3)+
  ylim(0,20)+
  xlab("Resistance") +
  ylab("Slope Mean") +
  ggtitle("Slope") +
  labs(tag="B")+
  theme_classic()+
  my_theme
slope


bat <- ggplot(data = mod_dat, aes(x = resistance, y = bat_mean)) +
  geom_boxplot(fill = "#D95F02", color = "black") +
  geom_jitter(width = 0.1, height = 0.3, alpha = 0.2, size=1) +
  ggsignif::geom_signif(comparisons = list(c("resistant", "transitioned")),
                        map_signif_level = TRUE,
                        tip_length = c(0.01, 0.01),
                        textsize=3)+
  ylim(5,23)+
  xlab("Resistance") +
  ylab("Depth (m) Mean") +
  ggtitle("Depth range") +
  labs(tag="")+
  theme_classic()+
  my_theme
bat


beuti <- ggplot(data = mod_dat, aes(x = resistance, y = beuti_month_obs)) +
  geom_boxplot(fill = "#7570B3", color = "black") +
  geom_jitter(width = 0.1, height = 0.3, alpha = 0.2, size=1) +
  ggsignif::geom_signif(comparisons = list(c("resistant", "transitioned")),
                        map_signif_level = TRUE,
                        y_position = 11,
                        tip_length = c(0.01, 0.01),
                        textsize=3)+
  ylim(2,13)+
  xlab("Resistance") +
  ylab("BEUTI Mean") +
  ggtitle("BEUTI") +
  labs(tag="")+
  theme_classic()+
  my_theme
beuti

sst <- ggplot(data = mod_dat, aes(x = resistance, y = sst_month_obs)) +
  geom_boxplot(fill = "#E7298A", color = "black") +
  geom_jitter(width = 0.1, height = 0.3, alpha = 0.2, size=1) +
  ggsignif::geom_signif(comparisons = list(c("resistant", "transitioned")),
                        map_signif_level = TRUE,
                        y_position = 15.3,
                        tip_length = c(0.01, 0.01),
                        textsize=3)+
  ylim(12,16)+
  xlab("Resistance") +
  ylab("SST Mean") +
  ggtitle("SST") +
  labs(tag="")+
  theme_classic()+
  my_theme
sst


kelp <- ggplot(data = mod_dat, aes(x = resistance, y = baseline_kelp)) +
  geom_boxplot(fill = "#66A61E", color = "black") +
  geom_jitter(width = 0.1, height = 0.3, alpha = 0.2, size=1) +
  ggsignif::geom_signif(comparisons = list(c("resistant", "transitioned")),
                        map_signif_level = TRUE,
                        tip_length = c(0.01, 0.01),
                        textsize=3)+
  ylim(50,220)+
  xlab("Resistance") +
  ylab("Kelp baseline density \n(no stipes per m²)") +
  ggtitle("Kelp baseline") +
  labs(tag="")+
  theme_classic()+
  my_theme
kelp

predictors <- ggpubr::ggarrange(slope, bat, kelp, beuti, sst, ncol=3, nrow=2) 
predictors


full_plot <- ggarrange(p1, predictors, nrow=1,  widths=c(1.5,2))

ggsave(full_plot, filename=file.path(figdir, "Fig5_predictors_new.png"), 
       width=7, height=6, bg="white", units="in", dpi=600)


