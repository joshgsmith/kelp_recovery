#Community data processing
#Joshua G. Smith
#April 25, 2023

rm(list=ls())
librarian::shelf(tidyverse, here)

################################################################################
#set directories and load data
basedir <- "/Volumes/seaotterdb$/kelp_recovery/"

kelp_fish_counts_raw <- read.csv(file.path(basedir, "data/subtidal_monitoring/raw/MLPA_kelpforest_fish.4.csv")) %>%
  janitor::clean_names()

kelp_upc_raw <- read.csv(file.path(basedir, "data/subtidal_monitoring/raw/MLPA_kelpforest_upc.4.csv")) %>%
  janitor::clean_names()
  
kelp_swath_raw <- read.csv(file.path(basedir, "data/subtidal_monitoring/raw/MLPA_kelpforest_swath.4.csv")) %>%
  janitor::clean_names()

kelp_taxon <- read.csv(file.path(basedir, "data/subtidal_monitoring/raw/MLPA_kelpforest_taxon_table.4.csv")) %>%
  janitor::clean_names() %>%
  dplyr::select(classcode, species_definition) %>%
  distinct()

site_table <- read.csv(file.path(basedir, "data/subtidal_monitoring/raw/MLPA_kelpforest_site_table.4.csv")) %>%
  janitor::clean_names() %>%
  dplyr::select(site, latitude, longitude, ca_mpa_name_short, mpa_class=site_designation, mpa_designation=site_status, 
                baseline_region)%>%
  distinct() #remove duplicates

################################################################################
#process kelp forest fish data

#select vars and clean
kelp_fish_build1 <- kelp_fish_counts_raw %>%
  dplyr::select(year, site, zone, level, transect, classcode, count)%>%
  group_by(year, site, zone, level, transect, classcode)%>%
  dplyr::summarize(total_count = sum(count)) #counts in raw data are grouped by size class. Take summary across all sizes

#join species names by class code 
kelp_fish_build2 <- left_join(kelp_fish_build1, kelp_taxon, by="classcode")

kelp_fish_build3 <- kelp_fish_build2 %>%
  dplyr::select(year, site, zone, level, transect, species=species_definition, total_count)%>%
  filter(!(species=='Sebastes atrovirens/carnatus/chrysomelas/caurinus'))


#add affiliated_mpa
kelp_fish_build4 <- left_join(kelp_fish_build3, site_table, by="site")

kelp_fish_build5 <- kelp_fish_build4 %>%
  ungroup() %>%
  dplyr::select(year, baseline_region, site,latitude, longitude, affiliated_mpa=ca_mpa_name_short,mpa_class, mpa_designation, zone, level, transect, species,
                total_count) %>%
  mutate(mpa_designation = ifelse(mpa_designation=="reference","ref",mpa_class)) %>%
  #make sure species are not duplicates by summarizing at the transect level (total counts)
  group_by(year, baseline_region, latitude, longitude, site, affiliated_mpa, mpa_class, mpa_designation, zone, level, transect, species) %>%
  dplyr::summarise(counts = sum(total_count))
  

#reshape to wide format 
kelp_fish_build6 <- kelp_fish_build5 %>%
  pivot_wider(names_from=species, values_from=counts, values_fill = 0) %>%
  janitor::clean_names()%>%
  #drop_na(region4)%>%
  mutate(MHW = ifelse(year>=2014 & year<=2016, "during",ifelse(year<2014, "before","after")))%>%
  dplyr::select(-c(no_organisms_present_in_this_sample)) %>%
  #filter central coast only
  filter(baseline_region == "CENTRAL") %>%
  #filter years >= 2007
  filter(year >= 2007)

unique(kelp_fish_build6$site)

#drop species and combine to higher taxa
#### see reference taxonomy table https://docs.google.com/spreadsheets/d/1vxy0XVOrlNhXD-i9tWL_F9G5h8S3S4OV/edit#gid=2031917236

kelp_fish_build7 <- as.data.frame(kelp_fish_build6) %>%
  
  #merge species
  dplyr::mutate(citharichthys_merge = rowSums(select(.,'citharichthys','citharichthys_sordidus','citharichthys_stigmaeus')
  ))%>%
  select(!(c('citharichthys','citharichthys_sordidus','citharichthys_stigmaeus'))) %>%
  #drop species
  select(!(c('anisotremus_davidsonii',
             'cheilotrema_saturnum',
             'unidentified_fish',
             'apodichthys_flavidus',
             'atherinopsidae',
             'aulorhynchus_flavidus',
             'bothidae',
             'clupeiformes',
             'cryptacanthodes_giganteus',
             'embiotocidae',
             'engraulis_mordax',
             'gobiesox_maeandricus',
             'hyperprosopon_anale',
             'hyperprosopon_argenteum',
             'hyperprosopon_ellipticum',
             'leiocottus_hirundo',
             'lethops_connectens',
             'neoclinus_blanchardi',
             'pholidae',
             'pleuronectidae',
             'rathbunella_alleni',
             'rathbunella_hypoplecta',
             'ronquilus_jordani',
             'sardinops_sagax',
             'sebastes',
             'stichaeidae',
             'syngnathus',
             'trachurus_symmetricus',
             'ulvicola_sanctaerosae')))%>%
  #drop string
  rename(citharichthys=citharichthys_merge) %>%
  dplyr::select(year, MHW, everything())


#Export
#write.csv(kelp_fish_build7,file.path(basedir, "/data/subtidal_monitoring/processed/kelp_fish_counts_CC.csv"), row.names = FALSE)


################################################################################
#process kelp swath

#select vars and clean
kelp_swath_build1 <- kelp_swath_raw %>%
  #drop species
  dplyr::filter(!(classcode == "STRPURREC" |
                    classcode == "MESFRAREC"))%>%
  dplyr::select(year, site, zone, transect, classcode, count, size)%>%
  group_by(year, site, zone, transect, classcode)%>%
  dplyr::summarize(total_count = sum(count), #counts in raw data are grouped by size class. Take summary across all sizes
                   total_size = sum(size)) %>% #this is for stipe counts only 
  mutate(total_count = ifelse(classcode == "MACPYRAD",total_size, total_count)) %>% #replace num plants with total stipes
  dplyr::select(!(total_size))
  
#join species names by class code 
kelp_swath_build2 <- left_join(kelp_swath_build1, kelp_taxon, by="classcode")

kelp_swath_build3 <- kelp_swath_build2 %>%
  dplyr::select(year, site, zone, transect, species=species_definition, total_count)


#add affiliated_mpa

kelp_swath_build4 <- left_join(kelp_swath_build3, site_table, by="site")

kelp_swath_build5 <- kelp_swath_build4 %>%
  ungroup() %>%
  dplyr::select(year, baseline_region, site,latitude, longitude, affiliated_mpa=ca_mpa_name_short,mpa_class, mpa_designation, zone, transect, species,
                total_count) %>%
  mutate(mpa_designation = ifelse(mpa_designation=="reference","ref",mpa_class)) %>%
  #make sure species are not duplicates by summarizing at the transect level (total counts)
  group_by(year, baseline_region, latitude, longitude, site, affiliated_mpa, mpa_class, mpa_designation, zone,  transect, species) %>%
  dplyr::summarise(counts = sum(total_count))



#reshape to wide format 
kelp_swath_build6 <- kelp_swath_build5 %>%
  pivot_wider(names_from=species, values_from=counts, values_fill = 0) %>%
  janitor::clean_names()%>%
  #drop_na(region4)%>%
  mutate(MHW = ifelse(year>=2014 & year<=2016, "during",ifelse(year<2014, "before","after")))%>%
  dplyr::select(-c(no_organisms_present_in_this_sample)) %>%
  #filter central coast only
  filter(baseline_region == "CENTRAL") %>%
  #filter years >= 2007
  filter(year >= 2007)

unique(kelp_swath_build6$site)


#drop species and combine to higher taxa
#### see reference taxonomy table https://docs.google.com/spreadsheets/d/1vxy0XVOrlNhXD-i9tWL_F9G5h8S3S4OV/edit#gid=2031917236

kelp_swath_build7 <- as.data.frame(kelp_swath_build6) %>%
  #merge species
  dplyr::mutate(urticina_merge = rowSums(dplyr::select(.,'urticina', 
                                                'urticina_coriacea', 'urticina_crassicornis', 'urticina_piscivora')
  ))%>%
  dplyr::select(!(c('urticina', 
             'urticina_coriacea', 'urticina_crassicornis', 'urticina_piscivora'))) %>%
  #drop species
  dplyr::select(!(c('anthopleura',
             'apostichopus',
             'asteroidea',
             'haliotis',
             'laminaria',
             'laminariales',
             'octopus',
             'pisaster',
             'pugettia_spp',
             'stylasterias_forreri',
             'urticina_columbiana',
             'urticina_columbiana_mcpeaki',
  )))%>%
  #drop string
  rename(urticina=urticina_merge)


#drop species that were never observed on the central coast
kelp_swath_zero_drop <- kelp_swath_build7 %>% filter(baseline_region=='CENTRAL') %>%
  dplyr::select(!(where(~ any(. != 0))))

kelp_swath_build8 <- kelp_swath_build7 %>% filter(baseline_region=='CENTRAL') %>%
  dplyr::select(where(~ any(. != 0))) %>%
  dplyr::select(year, MHW, everything())


names(kelp_swath_build8)

#Export
#write.csv(kelp_swath_build8,file.path(basedir, "/data/subtidal_monitoring/processed/kelp_swath_counts_CC.csv"), row.names = FALSE)


################################################################################
#process kelp upc

kelp_upc_build1 <- kelp_upc_raw%>%
  filter(category=="COVER")

#select vars and clean
kelp_upc_build2 <- kelp_upc_build1 %>%
  dplyr::select(year, site, zone, transect, classcode, count, pct_cov)%>%
  group_by(year, site, zone, transect, classcode)

#join species names by class code 
kelp_upc_build3 <- left_join(kelp_upc_build2, kelp_taxon, by="classcode")


kelp_upc_build4 <- kelp_upc_build3 %>%
  dplyr::select(year, site, zone, transect, classcode, species=species_definition, count, pct_cov)

#drop data 
kelp_upc_build5 <- kelp_upc_build4 %>%
  filter(!(species=="Haliotis"|
             species=="Laminaria"|
             species=="laminariales"|
             species=="UCSC 2000 Definition Of 'Other Red' Included All Red Algae Except Rhospp, Bushy, Leafy, Crucor And Erecor"|
             species=="UCSB 2000 Definition Of 'Other Red' Included All Red Algae Except Rhospp, Bushy, Leafy, Crucor And Erecor"|
             species=="Red Alga With Cylindrical Branches. Definition Used In 2000"|
             species=="UCSC and UCSB 1999 Definition Of 'Other Red' Included All Red Algae Except Rhospp, Crucor And Erecor"|
             species=="Laminariales Holdfast (Alive)"|
             species=="Sargassum"|
             species=="Aglaophenia struthionides"|
             species=="Macrosystis Holdfast (Dead)"|
             species=="Bare Rock"|
             species=="Bare Sand"|
             species=="Unidentified Fish"|
             species=="Actiniaria"|
             species=="Shell Debris"|
             species=="Sediment/Mud"|
             species=="Stylantheca papillosa"|
             species=="Diaperoforma californica"|
             species==""))


#recalculate percent cov
kelp_upc_build6 <- kelp_upc_build5 %>%
  group_by(year, site, zone, transect, classcode, species)%>%
  dplyr::summarize(sum_count = sum(count)) 

#check to make sure UPC counts add up to ~30 per transect               
check_counts <- kelp_upc_build5 %>%
  group_by(year, site, zone, transect)%>%
  dplyr::summarize(transect_total = sum(count))

#join percent cov with transect total
kelp_upc_build7 <- left_join(kelp_upc_build6, check_counts, by=c("year","site","zone","transect"))%>%
  mutate(pct_cov = (sum_count / transect_total)*100)%>%
  mutate_at(vars(pct_cov), round, 1)




#add affiliated_mpa

kelp_upc_build8 <- left_join(kelp_upc_build7, site_table, by="site")

kelp_upc_build9 <- kelp_upc_build8 %>%
  ungroup() %>%
  dplyr::select(year, baseline_region, site,latitude, longitude, affiliated_mpa=ca_mpa_name_short,mpa_class, mpa_designation, zone, transect, species,
                #sum_count, transect_total,
                pct_cov) %>%
  mutate(mpa_designation = ifelse(mpa_designation=="reference","ref",mpa_class)) 



#reshape to wide format
kelp_upc_build10 <- kelp_upc_build9 %>%
  pivot_wider(names_from=species, values_from=pct_cov) %>%
  janitor::clean_names()%>%
  #drop_na(region4) %>%
  mutate(MHW = ifelse(year>=2014 & year<=2016, "during",ifelse(year<2014, "before","after")))%>%
  #dplyr::select(year, region3, region4, affiliated_mpa, mpa_defacto_class, MHW, everything())%>%
  #filter central coast only
  filter(baseline_region == "CENTRAL") %>%
  #filter years >= 2007
  filter(year >= 2007) %>%
  mutate_at(c(11:67), ~replace_na(.,0)) %>%
  dplyr::select(year, MHW, everything())
  #filter(region3 == 'central') %>%
  #select(where(~ any(. != 0)))
  #mutate_all(~ifelse(is.nan(.), NA, .))

#Export
#write.csv(kelp_upc_build10,file.path(basedir, "/data/subtidal_monitoring/processed/kelp_upc_cov_CC.csv"), row.names = FALSE)




