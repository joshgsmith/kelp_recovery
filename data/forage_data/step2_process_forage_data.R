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
         kelptype = tolower(kelptype),
         success = tolower(success),
         sb_area = toupper(sb_area)
         
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

