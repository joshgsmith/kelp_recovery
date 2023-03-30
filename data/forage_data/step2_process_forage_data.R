#Joshua G. Smith
#jossmith@mbayaq.org
#Script initiated March 29, 2023

rm(list=ls())

#load packages
librarian::shelf(tidyverse, readxl, here, DataExplorer)

#set directories 
datin <- "/Volumes/seaotterdb$/kelp_recovery/data/foraging_data/raw"
datout <- "/Volumes/seaotterdb$/kelp_recovery/data/foraging_data/processed"
figdir <- here::here("analyses","figures")

#read 2016-current data

for_dat <- read_xlsx(file.path(datin, "Forage_data_2016tocurrent.xlsx"))
dives <- read_xlsx(file.path(datin, "Forage_dives_2016tocurrent.xlsx"))
index <- read_xlsx(file.path(datin, "Forage_indexi_2016tocurrent.xlsx"))

################################################################################
#explore data

str(dat)
str(dives)
str(index)

plot_intro(dat)
plot_intro(dives)
plot_intro(index)

################################################################################
#process foraging data

for_dat_build1 <- for_dat %>% janitor::clean_names() %>%
                    #correct inconsistencies
                    mutate(qualifier = factor(tolower(qualifier)),
                           pup_behav = factor(tolower(pup_behav)),
                           mom_resp = factor(tolower(mom_resp)),
                           outcome = factor(tolower(outcome)),
                           prey = factor(prey),
                           tooltype = factor(tooltype),
                           pupsh = factor(pupsh),
                           pup_behav = factor(pup_behav),
                           stolenfr = factor(stolenfr),
                           steal = factor(steal),
                           sizecm = factor(sizecm)
                          )

str(for_dat_build1)

#check for inconsistencies
unique(for_dat_build1$preynum)       
unique(for_dat_build1$prey)    
unique(for_dat_build1$number)    
unique(for_dat_build1$size)    
unique(for_dat_build1$qualifier) #make lower
unique(for_dat_build1$ht)
unique(for_dat_build1$tooltype)
unique(for_dat_build1$tool_number)
unique(for_dat_build1$pupsh)
unique(for_dat_build1$pupshare_number)
unique(for_dat_build1$pup_behav) #make lower
unique(for_dat_build1$mom_resp) #make lower
unique(for_dat_build1$outcome) #make lower
unique(for_dat_build1$stolenfr) 
unique(for_dat_build1$stolenfr_number)
unique(for_dat_build1$steal)
unique(for_dat_build1$steal_number)
unique(for_dat_build1$sizecm)
unique(for_dat_build1$gulls)


################################################################################
#process dives

dive_build1 <- dives %>% janitor::clean_names() %>%
  #correct inconsistencies
  mutate(where = tolower(where),
         #fix incorrect sign
         long_obs_deg = ifelse(long_obs_deg == 121,-121,long_obs_deg),
         #fix kelp type
         kelptype = ifelse(kelptype == "xx","x",kelptype),
         kelptype = factor(tolower(kelptype)),
         success = factor(tolower(success)),
         sb_area = factor(toupper(sb_area))
         
  ) 

colnames(dive_build1)

str(dive_build1)

#check for inconsistencies
unique(dive_build1$subbout)       
unique(dive_build1$where)  
unique(dive_build1$lat_obs_deg)  
unique(dive_build1$long_obs_deg)
unique(dive_build1$long_obs_min)
unique(dive_build1$lat_obs)
unique(dive_build1$long_obs)
unique(dive_build1$long)
unique(dive_build1$shore)
unique(dive_build1$depth)
unique(dive_build1$canopy)
unique(dive_build1$kelptype)
unique(dive_build1$divenum)
unique(dive_build1$dt)
unique(dive_build1$st)
unique(dive_build1$success)
unique(dive_build1$sb_area)



################################################################################
#process index

index_build1 <- index %>% janitor::clean_names() %>%
                mutate(visib = factor(tolower(visib)),
                       sex = factor(sex),
                       ageclass = factor(ageclass),
                       status = factor(status),
                       pupsz = factor(pupsz),
                       consort = factor(consort),
                       obsbeg = factor(obsbeg),
                       obsend = factor(obsend),
                       daynight = factor(daynight),
                       area = factor(area),
                       sky = factor(sky),
                       winddir = factor(winddir),
                       seaopen = factor(seaopen),
                       seafix = factor(seafix),
                       swell = factor(swell)
                       )

str(index_build1)

colnames(index_build1)

unique(index_build1$sex)  
unique(index_build1$ageclass)  
unique(index_build1$status)  
unique(index_build1$pupsz)  
unique(index_build1$ageweeks)  
unique(index_build1$agequal)  
unique(index_build1$consort)  
unique(index_build1$gen_locat)  
unique(index_build1$atos)  
unique(index_build1$area)  
unique(index_build1$fixtype) 
unique(index_build1$sky)  
unique(index_build1$rain)  
unique(index_build1$windspeed)  
unique(index_build1$winddir)  
unique(index_build1$temperature)  
unique(index_build1$seaopen)  
unique(index_build1$seafix)  
unique(index_build1$swell)  
unique(index_build1$visib)  



################################################################################
#merge data

forage_join <- for_dat_build1 %>%
              #filter(steal == "No")%>% 
                #select variables of interest
              dplyr::select(foragdiv_id, foragdata_id, preynum, prey, number,
                            size, qualifier)%>%
              mutate_if(is.character, str_trim)

dive_join <- dive_build1 %>%
             dplyr::select(foragdiv_id, bout, subbout, lat, long, canopy,
                           kelptype, divenum, dt, success) %>%
              mutate_if(is.character, str_trim) %>%
              #filter successful dives only
              filter(success == "y")

index_join <- index_build1 %>%
              dplyr::select(bout, date, otterno, sex, ageclass, status)%>%
             mutate_if(is.character, str_trim)


#merge forage data and bout info
data_build1 <- left_join(dive_join, forage_join, by="foragdiv_id")
anti_build <- anti_join(dive_join, forage_join, by="foragdiv_id")


#merge with index
data_build2 <- left_join(data_build1, index_join, by="bout")
anti_build2 <- anti_join(data_build1, index_join, by="bout")


################################################################################
#process data build

data_build3 <- data_build2 %>%
                mutate(year = as.numeric(format(date,'%Y')),
                       month = as.numeric(format(date,'%m')),
                                          day = as.numeric(format(date,'%d')))%>%
                dplyr::select(year, month, day, foragdiv_id, foragdata_id, bout, subbout,
                              lat, long, otterno, everything()
                              ) %>%
                filter(year >= 2016)


write.csv(data_build3, file.path(datout, "foraging_data_2016_2023.csv"), row.names=FALSE)




