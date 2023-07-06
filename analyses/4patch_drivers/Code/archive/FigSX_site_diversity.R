#Community data processing
#Joshua G. Smith
#April 28, 2023

rm(list=ls())
librarian::shelf(tidyverse, here, vegan, cowplot, ggpubr)

library(dplyr)

################################################################################
#set directories and load data
basedir <- "/Volumes/seaotterdb$/kelp_recovery/"
figdir <- here::here("analyses","4patch_drivers","Figures")

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
#process swath data

swath_dat <- swath_raw %>% dplyr::select(12:ncol(.))
swath_groups <- swath_raw %>% dplyr::select(1:11)

swath_richness <- data.frame(S.obs = apply(swath_dat[,1:59]>0, 1, sum))
swath_evenness <- diversity(swath_dat)/log(specnumber(swath_dat))
swath_shannon <- diversity(swath_dat, index="shannon")
swath_simpson <- diversity(swath_dat, index="simpson")
swath_abund <- rowSums(swath_dat[,1:59])

swath_alphadiv <- cbind(swath_groups, swath_richness, swath_shannon, swath_simpson, swath_evenness, swath_abund)%>%
  mutate(MHW = str_to_sentence(MHW),
         MHW = factor(MHW, levels=c("Before","During","After")))


#process upc data

upc_dat <- upc_raw %>% dplyr::select(12:ncol(.))
upc_groups <- upc_raw %>% dplyr::select(1:11)

total_area_covered <- rowSums(upc_dat)

abundance <- upc_dat / total_area_covered * 166.67 #convert %cov to abundance, multiply by scalar, in this case 60m2 transect, 166.67 (which is 10,000 divided by 60)

upc_dat[,1:56] <- abundance

#convert perc cov to relative abundance

upc_richness <- data.frame(S.obs = apply(upc_dat[,1:56]>0, 1, sum))
upc_evenness <- diversity(upc_dat)/log(specnumber(upc_dat))
upc_shannon <- diversity(upc_dat, index="shannon")
upc_simpson <- diversity(upc_dat, index="simpson")
upc_abund <- rowSums(upc_dat[,1:56])

upc_alphadiv <- cbind(upc_groups, upc_richness, upc_shannon, upc_simpson, upc_evenness, upc_abund)%>%
  mutate(MHW = str_to_sentence(MHW),
         MHW = factor(MHW, levels=c("Before","During","After")))


#process fish data

fish_dat <- fish_raw %>% dplyr::select(13:ncol(.))
fish_groups <- fish_raw %>% dplyr::select(1:12)

fish_richness <- data.frame(S.obs = apply(fish_dat[,1:111]>0, 1, sum))
fish_evenness <- diversity(fish_dat)/log(specnumber(fish_dat))
fish_shannon <- diversity(fish_dat, index="shannon")
fish_simpson <- diversity(fish_dat, index="simpson")
fish_abund <- rowSums(fish_dat[,1:111])

fish_alphadiv <- cbind(fish_groups, fish_richness, fish_shannon, fish_simpson, fish_evenness, fish_abund)%>%
  mutate(MHW = str_to_sentence(MHW),
         MHW = factor(MHW, levels=c("Before","During","After")),
         fish_evenness = ifelse(fish_evenness == "NaN",NA,fish_evenness))
  
################################################################################
###test for significant differences for swath


# Load the data
my_data <- swath_alphadiv

# Get a vector of unique site names
site_names <- unique(my_data$site)


# Create an empty data frame to store the results
results <- data.frame(site = character(),
                      intercept = numeric(),
                      slope = numeric(),
                      p_value = numeric(),
                      stringsAsFactors = FALSE)

# Loop over each site and fit a linear regression model
for (site_name in site_names) {
  site_data <- subset(my_data, site == site_name)
  model <- lm(swath_shannon ~ year, data = site_data)
  p_value <- summary(model)$coefficients[2,4]
  p_value <- format(p_value, scientific = FALSE)
  results <- rbind(results, data.frame(site = site_name,
                                       intercept = coef(model)[1],
                                       slope = round(coef(model)[2],2),
                                       p_value = round(as.numeric(p_value),3)))
}

View(results)

swath_sig <- results %>%
                   mutate(p_new = ifelse(p_value == 0 ,0.001,p_value),
                     sig = ifelse(p_new < 0.05,"*",""),
                   p_lab = paste0("p < ",p_new,sig),
                   slope = paste("Linear slope:",slope)) %>% distinct()

################################################################################
###test for significant differences for fish
                  

# Load the data
my_data <- fish_alphadiv

# Get a vector of unique site names
site_names <- unique(my_data$site)


# Create an empty data frame to store the results
results <- data.frame(site = character(),
                      intercept = numeric(),
                      slope = numeric(),
                      p_value = numeric(),
                      stringsAsFactors = FALSE)

# Loop over each site and fit a linear regression model
for (site_name in site_names) {
  site_data <- subset(my_data, site == site_name)
  model <- lm(fish_shannon ~ year, data = site_data)
  p_value <- summary(model)$coefficients[2,4]
  p_value <- format(p_value, scientific = FALSE)
  results <- rbind(results, data.frame(site = site_name,
                                       intercept = coef(model)[1],
                                       slope = round(coef(model)[2],2),
                                       p_value = round(as.numeric(p_value),3)))
}

View(results)

fish_sig <- results %>%
  mutate(p_new = ifelse(p_value == 0 ,0.001,p_value),
         sig = ifelse(p_new < 0.05,"*",""),
         p_lab = paste0("p < ",p_new,sig),
         slope = paste("Linear slope:",slope)) %>% distinct()

################################################################################
###test for significant differences for UPC


# Load the data
my_data <- upc_alphadiv

# Get a vector of unique site names
site_names <- unique(my_data$site)


# Create an empty data frame to store the results
results <- data.frame(site = character(),
                      intercept = numeric(),
                      slope = numeric(),
                      p_value = numeric(),
                      stringsAsFactors = FALSE)

# Loop over each site and fit a linear regression model
for (site_name in site_names) {
  site_data <- subset(my_data, site == site_name)
  model <- lm(upc_shannon ~ year, data = site_data)
  p_value <- summary(model)$coefficients[2,4]
  p_value <- format(p_value, scientific = FALSE)
  results <- rbind(results, data.frame(site = site_name,
                                       intercept = coef(model)[1],
                                       slope = round(coef(model)[2],2),
                                       p_value = round(as.numeric(p_value),3)))
}

View(results)

upc_sig <- results %>%
  mutate(p_new = ifelse(p_value == 0 ,0.001,p_value),
         sig = ifelse(p_new < 0.05,"*",""),
         p_lab = paste0("p < ",p_new,sig),
         slope = paste("Linear slope:",slope)) %>% distinct()

################################################################################
###test for significant differences overall

#swath
  model <- lm(swath_shannon ~ year, data = swath_alphadiv)
  p_value <- summary(model)$coefficients[2,4]
  p_value <- format(p_value, scientific = FALSE)
  swath_results <- cbind(data.frame(intercept = coef(model)[1],
                                       slope = round(coef(model)[2],2),
                                       p_value = round(as.numeric(p_value),3))) %>%
        mutate(p_new = ifelse(p_value == 0 ,0.001,p_value),
         sig = ifelse(p_new < 0.05,"*",""),
         p_lab = paste0("p < ",p_new,sig),
         slope = paste("Linear slope:",slope)) %>% distinct()
  
#Fish
  model <- lm(fish_shannon ~ year, data = fish_alphadiv)
  p_value <- summary(model)$coefficients[2,4]
  p_value <- format(p_value, scientific = FALSE)
  fish_results <- cbind(data.frame(intercept = coef(model)[1],
                                    slope = round(coef(model)[2],2),
                                    p_value = round(as.numeric(p_value),3))) %>%
    mutate(p_new = ifelse(p_value == 0 ,0.001,p_value),
           sig = ifelse(p_new < 0.05,"*",""),
           p_lab = paste0("p < ",p_new,sig),
           slope = paste("Linear slope:",slope)) %>% distinct()

  #UPC
  model <- lm(upc_shannon ~ year, data = upc_alphadiv)
  p_value <- summary(model)$coefficients[2,4]
  p_value <- format(p_value, scientific = FALSE)
  upc_results <- cbind(data.frame(intercept = coef(model)[1],
                                   slope = round(coef(model)[2],2),
                                   p_value = round(as.numeric(p_value),3))) %>%
    mutate(p_new = ifelse(p_value == 0 ,0.001,p_value),
           sig = ifelse(p_new < 0.05,"*",""),
           p_lab = paste0("p < ",p_new,sig),
           slope = paste("Linear slope:",slope)) %>% distinct()

################################################################################
#plot site level

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


all_three <- ggplot()+
  #swath
  geom_point(aes(x = year, y=swath_shannon, color = "Swath"), alpha = 0.2, size=0.5, data = swath_alphadiv %>% mutate(site = gsub("_", " ", site)),
             position = position_jitter(width = 0.02)) +
  geom_smooth(method = "auto", se = TRUE, size = 0.5, alpha = 0.3, aes(x = year, y=swath_shannon, color = "Swath"),data = swath_alphadiv %>% mutate(site = gsub("_", " ", site))) +
  geom_text(aes(x=2018, y=2.8, label=p_lab, color="Swath"), size=3, data=swath_sig %>% mutate(site = gsub("_", " ", site)))+ #add p
  geom_text(aes(x=2014, y=2.8, label=slope, color="Swath"), size=3, data=swath_sig %>% mutate(site = gsub("_", " ", site)))+ #add slope
  #fish
  geom_point(aes(x = year-0.1, y=fish_shannon, color = "Fish"), alpha = 0.2, size=0.5, data = fish_alphadiv %>% mutate(site = gsub("_", " ", site)), 
             position = position_jitter(width = 0.02)) +
  geom_smooth(method = "auto", se = TRUE, size = 0.5, alpha = 0.3, aes(x = year, y=fish_shannon, color = "Fish"),data = fish_alphadiv %>% mutate(site = gsub("_", " ", site)) ) +
  geom_text(aes(x=2018, y=2.6, label=p_lab, color="Fish"), size=3, data=fish_sig %>% mutate(site = gsub("_", " ", site)))+ #add p
  geom_text(aes(x=2014, y=2.6, label=slope, color="Fish"), size=3, data=fish_sig %>% mutate(site = gsub("_", " ", site)))+ #add slope
  #upc
  geom_point(aes(x = year+0.1, y=upc_shannon, color = "UPC"), alpha = 0.2, size=0.5, data = upc_alphadiv %>% mutate(site = gsub("_", " ", site)), 
             position = position_jitter(width = 0.02)) +
  geom_smooth(method = "auto", se = TRUE, size = 0.5, alpha = 0.3, aes(x = year, y=upc_shannon, color = "UPC"),data = upc_alphadiv %>% mutate(site = gsub("_", " ", site))) +
  geom_text(aes(x=2018, y=2.4, label=p_lab, color="UPC"), size=3, data=upc_sig %>% mutate(site = gsub("_", " ", site)))+ #add p
  geom_text(aes(x=2014, y=2.4, label=slope, color="UPC"), size=3, data=upc_sig %>% mutate(site = gsub("_", " ", site)))+ #add slope
  facet_wrap(~site)+
  #color
  scale_color_manual(values = c("#009E73", "#E69F00","#CC79A7"), labels = c("Fish", "Swath","UPC"))+ 
  theme_bw() + my_theme

all_three



#ggsave(combined_plot, filename=file.path(figdir, "Fig2_diversity_overall.png"), bg = "white",
 #      width=6.5, height=6, units="in", dpi=600)

























