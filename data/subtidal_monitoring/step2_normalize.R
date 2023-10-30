#Community data processing
#Joshua G. Smith
#April 28, 2023

rm(list=ls())
librarian::shelf(tidyverse, here, vegan)


################################################################################
#set directories and load data
basedir <- "/Volumes/seaotterdb$/kelp_recovery/"

fish_raw <- read.csv(file.path(basedir, "data/subtidal_monitoring/processed/kelp_fish_counts_CC.csv")) 

upc_raw <- read.csv(file.path(basedir, "data/subtidal_monitoring/processed/kelp_upc_cov_CC.csv")) 

swath_raw <- read.csv(file.path(basedir, "data/subtidal_monitoring/processed/kelp_swath_counts_CC.csv")) 


################################################################################
#pre-processing

#drop species that were never encountered
fish_build1 <- fish_raw %>% dplyr::select(where(~ any(. != 0)))
upc_build1 <- upc_raw %>% dplyr::select(where(~ any(. != 0)))
swath_build1 <- swath_raw %>% dplyr::select(where(~ any(. != 0)))

#check the number of transects per site. We are going to take the mean
t_fish <- fish_build1 %>% group_by(year, site) %>% summarize(n_t = n())
ggplot(data = t_fish, aes(x = year, y=n_t)) +
  geom_bar(stat = "identity") +
  facet_wrap(~site) +
  labs(x = "Year", y = "Count") +
  ggtitle("Bar plot of Count by Year and Site")


#take mean spp counts per transect
fish_build2 <- fish_build1 %>%
  dplyr::group_by(year, MHW, baseline_region, latitude, longitude, site,
                  affiliated_mpa, mpa_class, mpa_designation) %>%
  dplyr::summarize(across(4:55, mean, na.rm = TRUE))


upc_build2 <- upc_build1 %>%
  dplyr::group_by(year, MHW, baseline_region, latitude, longitude, site,
                  affiliated_mpa, mpa_class, mpa_designation) %>%
  #remove UPC species that are recorded by swath
  dplyr::summarize(across(3:41, mean, na.rm = TRUE)) %>%
  dplyr::select(!c('macrocystis_pyrifera',
                   'stylaster_californicus',
                   'stephanocystis_osmundacea')) 


swath_build2 <- swath_build1 %>%
  dplyr::group_by(year, MHW, baseline_region, latitude, longitude, site,
                  affiliated_mpa, mpa_class, mpa_designation) %>%
  dplyr::summarize(across(3:59, mean, na.rm = TRUE))


################################################################################
#standardize 

##kelp fish
fish_stan <- decostand(fish_build2[10:61], method = "max",margin = 2, na.rm=TRUE)
nrow(fish_build2) #check length
nrow(fish_stan) #check length
#join with group vars
fish_stan1 <- cbind(fish_build2[1:9],fish_stan)

###kelp upc
upc_stan <- decostand(upc_build2[10:45], method = "max",margin = 2, na.rm=TRUE)
nrow(upc_build2) #check length
nrow(upc_stan) #check length
#join with group vars
upc_stan1 <- cbind(upc_build2[1:9],upc_stan) 
  

###kelp swath
swath_stan <- decostand(swath_build2[10:66], method = "max",margin = 2, na.rm=TRUE)
nrow(swath_build2) #check length
nrow(swath_stan) #check length

#join with group vars
swath_stan1 <- cbind(swath_build2[1:9],swath_stan)



################################################################################
#merge

kelp_stan_all <- swath_stan1 %>%
  ###NOTE: joning based on swath_stan, which has the data on urchins and kelp
  left_join(upc_stan1, by = c("year", "MHW", "baseline_region", "latitude", "longitude", "site", "affiliated_mpa", "mpa_class", "mpa_designation")) %>%
  left_join(fish_stan1, by = c("year", "MHW", "baseline_region", "latitude", "longitude", "site", "affiliated_mpa", "mpa_class", "mpa_designation"))

#check if it worked
nrow(swath_stan1)
nrow(kelp_stan_all)

#export
write.csv(kelp_stan_all,file.path(here::here("analyses","4patch_drivers","Output"),"kelp_stan_CC.csv"), row.names = FALSE)

#last write 30 Oct 2023





