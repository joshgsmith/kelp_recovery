#Communtiy change analyses
#Joshua G. Smith; jossmith@mbayaq.org

rm(list=ls())

librarian::shelf(tidyverse, here, vegan, ggplot2, cluster, ggforce)


################################################################################
#set directories and load data
basedir <- "/Volumes/seaotterdb$/kelp_recovery/"
figdir <- here::here("analyses","4patch_drivers","Figures")

#load standardized dat
stan_dat <- read.csv(file.path(basedir, "data/subtidal_monitoring/processed/kelp_stan_CC.csv")) 

################################################################################
#process standardized dat

#replace any NAs with 0
stan_dat <- stan_dat %>% mutate(across(where(is.numeric), ~replace_na(., 0))) %>%
  #select sites in Carmel and Monterey Bay only
  dplyr::filter(latitude >= 36.46575 & latitude <= 36.64045) %>%
  #drop sea stars and purple urchins
  #dplyr::select(!(c(strongylocentrotus_purpuratus,
  #                 pycnopodia_helianthoides,
  #                macrocystis_pyrifera)))%>%
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
#Define grouping variables for data matrix

#define group vars
stan_group_vars <- stan_dat %>%
                    #create identifier
                    mutate(identifier = paste(site, year)) %>%
                    dplyr::select(identifier, everything())

#define data for matrix
stan_mat_data <- stan_group_vars %>% ungroup() %>% dplyr::select(identifier, 11:ncol(.))

#convert identifier to row name
rownames(stan_mat_data) <- stan_mat_data[,1]
stan_mat_data <- stan_mat_data[,-1]

stan_identifier <- stan_group_vars %>% dplyr::select(identifier)


################################################################################
#Create dissimilarity matrix


bc_mat <- as.matrix(vegan::vegdist(stan_mat_data, method = "bray"))

#join group vars
dist_dat <- cbind(stan_identifier, bc_mat)

#create header names to match square matrix
colnames(dist_dat)[2:ncol(dist_dat)] <- dist_dat[,1]

#convert to three column format
dist_long <- setNames(melt(dist_dat),c('site_year1','site_year2','dissim'))


################################################################################
#Format matrix

dist_long1 <- dist_long %>%
          data.frame()%>%
          #bust apart the identifier 
          mutate(site_1 = as.factor(as.character(word(site_year1,1,-2))),
                 site_2 = as.factor(as.character(word(site_year2,1,-2))),
                 year_1 = as.factor(as.character(word(site_year1, -1))),
                 year_2 = as.factor(as.character(word(site_year2, -1)))) %>%
          #remove white space
          mutate(site_1 = trimws(as.factor(str_trim(site_1))),
            site_2 = trimws(as.factor(str_trim(site_2)))) %>%
          #filter by year = year for between site dissim
          filter(year_1 == year_2) %>%
          #rearrange and drop identifier
          dplyr::select(site_1, site_2, year_1, year_2, dissim) %>%
          group_by(year_1, year_2, site_1, site_2, dissim)%>%
          #drop within site comparisons
          dplyr::filter(!(site_1 == site_2)) %>%
          #take lower triangle  but removing duplicates
          distinct(dissim)
  

################################################################################
#aggregate data

dist_long2 <- dist_long1 %>%
              group_by(year_1, site_1)%>%
              dplyr::summarize(n_dissim = n(),
                               mean_dissim = mean(dissim, na.rm=TRUE),
                               sd_dissim = sd(dissim, na.rm=TRUE)
                               )











