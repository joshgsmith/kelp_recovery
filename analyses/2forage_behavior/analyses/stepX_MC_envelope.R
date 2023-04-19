#Joshua G. Smith, jossmith@mbayaq.org

rm(list=ls())

#required packages
librarian::shelf(tidyverse, maptools, rgdal, spatstat, sf)

#set directories 
basedir <- "/Volumes/seaotterdb$/kelp_recovery/data"
figdir <- here::here("analyses","figures")


#read foraging data
forage_dat <- read.csv(file.path(basedir, "foraging_data/processed/foraging_data_2016_2023.csv")) 

#read forage metadata
forage_meta <- read.csv(file.path(basedir, "foraging_data/processed/forage_metadata.csv")) 


################################################################################
#Step 1 -- process foraging dat for join

forage_build1 <- forage_dat %>%
  mutate(quarter = ifelse(month == 1 | month == 2 | month == 3, 1,
                          ifelse(month == 4 | month == 5 | month == 6,2,
                                 ifelse(month == 7 | month == 8 | month == 9,3,4)))) %>%
  #filter urchin prey only
  filter(prey == "pur") 

focal_patch <-  forage_dat %>%
  mutate(quarter = ifelse(month == 1:3,1,
                          ifelse(month == 4:6,2,
                                 ifelse(month==7:9,3,4)))) %>%
  #filter urchin prey only
  filter(prey == "pur") %>%
  #define focal patch
  group_by(bout) %>%
  dplyr::summarize(n_dives = n()) %>%
  dplyr::mutate(focal_patch = ifelse(n_dives > 3, "yes","no"))


#join

forage_build2 <- left_join(forage_build1, focal_patch, by="bout")



#aggregate data by bout

forage_build3 <- forage_build2 %>%
  group_by(year, quarter, bout, focal_patch) %>%
  summarize(n_prey = sum(preynum),
            lat = mean(lat),
            long = mean(long)
  )


################################################################################
#Step 2 -- build mante carlo


focal_patch <- forage_build3 %>% filter(focal_patch == "yes") %>% ungroup() %>%
                  #drop NA
                  filter(!(is.na(lat)))
  


focal_sp <- SpatialPointsDataFrame(coords = focal_patch[,c("long", "lat")], 
                                   data = focal_patch[,c("year", "quarter","bout","n_prey")], 
                                   proj4string = CRS("+init=epsg:3310"))


SP <- as(focal_sp, "SpatialPoints")
P <- as(SP, "ppp")

E <- envelope(P, Kest,nsim = 500)
plot(E, main = "MCMC Global envelopes \nfocal patch (>3 successful dives)")













