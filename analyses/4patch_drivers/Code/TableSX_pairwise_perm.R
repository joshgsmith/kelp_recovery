#Communtiy change analyses
#Joshua G. Smith; jossmith@mbayaq.org

rm(list=ls())

#install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")

librarian::shelf(tidyverse, here, vegan, ggplot2, cluster, ggforce, 
                 pairwiseAdonis, broom)




################################################################################
#set directories and load data
basedir <- "/Volumes/seaotterdb$/kelp_recovery/"

stan_dat <- read.csv(file.path(basedir, "data/subtidal_monitoring/processed/kelp_stan_CC.csv")) 

tab_dir <- here::here("analyses","4patch_drivers","Tables")

################################################################################
#pre-analysis processing

#replace any NAs with 0
stan_dat <- stan_dat %>% mutate(across(where(is.numeric), ~replace_na(., 0))) %>%
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
                    site == "BIRD_ROCK")) # %>%
  #test permova if urchins, kelp, and sea stars are removed
  #dplyr::select(!(c(strongylocentrotus_purpuratus, macrocystis_pyrifera,
   #                 pycnopodia_helianthoides, pisaster_brevispinus,
    #                pisaster_giganteus, pisaster_ochraceus,
     #              dermasterias_imbricata)))


################################################################################
#prepare data 

#----------------process standardized data--------------------------------------

#define group vars
stan_group_vars <- stan_dat %>% 
                   dplyr::select(1:9) %>%
                    mutate(site_period = paste(site, MHW),
                           outbreak_period = ifelse(year <2014, "before","after"),
                           site_outbreak = paste(site, ifelse(year <2014, "before","after")))

#define data for ordination
stan_ord_dat <- stan_dat %>% dplyr::select(10:ncol(.))

#standardize to max 
stan_rel <- decostand(stan_ord_dat, method = "hellinger")

#generate a BC mat with stan dat
stan_max_distmat <- vegdist(stan_rel, method = "bray", na.rm = T)


################################################################################
#build pairwise PERMANOVA for regional comparison

#pariwise permanova for regional
set.seed(1985)

pair_perm <- pairwise.adonis2(stan_max_distmat ~ outbreak_period, 
                                    data = stan_group_vars, 
                                    permutations = 999,
                                    num_cores = 20)

##### EXTRACT OUTPUT
# Remove the different list item
pair_perm$parent_call <- NULL

# make it a data frame
regional_result_tab <- pair_perm %>% map_dfr(~tidy(.x), .id="name")%>%
  filter(term == "outbreak_period")


################################################################################
#build pairwise PERMANOVA for site-level comparison

#pariwise permanova for regional
set.seed(1985)

pair_perm <- pairwise.adonis2(stan_max_distmat ~ site_outbreak, 
                              data = stan_group_vars, 
                              permutations = 999,
                              num_cores = 20)

##### EXTRACT OUTPUT
# Remove the different list item
pair_perm$parent_call <- NULL

# make it a data frame
site_result_tab <- pair_perm %>% map_dfr(~tidy(.x), .id="name")%>%
  filter(term == "site_outbreak") %>%
  mutate(group_1 = str_extract(name, ".*(?=_vs)"), #extract everything before "_vs"
         group_2 = str_extract(name, "(?<=vs_).*"), #extract everything after "vs_"
         site_1 = str_extract(group_1, "^[^ ]+"), #extract everything before " "
         site_2 = str_extract(group_2, "^[^ ]+"), #extract everything before " "
         period_1 = str_extract(group_1, "\\w+$"),
         period_2 = str_extract(group_2, "\\w+$")) %>%
  #filter matrix
  filter(site_1 == site_2,
         period_1 == "before")%>%
  dplyr::select(site = site_1,
                period_1, period_2, df, SumOfSqs, R2, statistic, p.value)



write.csv(site_result_tab, file.path(tab_dir, "TableS2_pairwise_permanova.csv"), row.names = FALSE)


