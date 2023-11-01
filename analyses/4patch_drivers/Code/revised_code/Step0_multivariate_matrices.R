#Communtiy change analyses
#Joshua G. Smith; jossmith@mbayaq.org

#create distance matrices and ordinate data

rm(list=ls())

librarian::shelf(tidyverse, here, vegan, ggplot2, cluster, ggforce)


################################################################################
#set directories and load data

figdir <- here::here("analyses","4patch_drivers","Figures")
basedir <- here::here("analyses","4patch_drivers","Output")

#load standardized dat
stan_dat <- read.csv(file.path(basedir, "kelp_stan_CC.csv")) 

#load fish
fish_raw <- read.csv(file.path(basedir, "kelp_fish_counts_CC.csv"))

#load upc
upc_raw <- read.csv(file.path(basedir, "kelp_upc_cov_CC.csv")) 

#load swath
swath_raw <- read.csv(file.path(basedir, "kelp_swath_counts_CC.csv")) 


################################################################################
#process data

#replace any NAs with 0
stan_dat <- stan_dat %>% mutate(across(where(is.numeric), ~replace_na(., 0))) 

#take the mean across replicate transects for a given site
fish_sum <- fish_raw %>% group_by(year, MHW, site) %>%
  dplyr::summarize(across(10:114, mean, na.rm = TRUE))

swath_sum <- swath_raw %>% group_by(year, MHW, site) %>%
  dplyr::summarize(across(9:67, mean, na.rm = TRUE))

upc_sum <- upc_raw %>% group_by(year, MHW, site) %>%
  dplyr::summarize(across(9:64, mean, na.rm = TRUE))


################################################################################
#prepare data for ordination

#----------------process standardized data--------------------------------------

#define group vars
stan_group_vars <- stan_dat %>% dplyr::select(1:9)

#define data for ordination
stan_ord_dat <- stan_dat %>% dplyr::select(10:ncol(.))

#standardize to max 
stan_rel <- decostand(stan_ord_dat, method = "hellinger")

#generate a BC mat with stan dat
stan_max_distmat <- vegdist(stan_rel, method = "bray", na.rm = T)

#generate a BC mat with ord dat
stan_untransformed_distmat <- vegdist(stan_ord_dat, method = "bray", na.rm = T)


################################################################################
#ordinate data

set.seed(1985)
num_cores = 8

#ordinate stan dat
stan_ord <- metaMDS(stan_max_distmat, distance = "bray", parallel = num_cores, trymax=300)
stan_untrans_ord <- metaMDS(stan_untransformed_distmat, distance = "bray", parallel = num_cores, trymax=300)


save(file = paste(file.path(basedir,"multivariate_data.Rdata")),stan_dat, stan_group_vars, stan_ord,
   stan_ord_dat, stan_rel, stan_max_distmat)

#last write 31 Oct 2023



