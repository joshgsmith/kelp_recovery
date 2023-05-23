#Communtiy change analyses
#Joshua G. Smith; jossmith@mbayaq.org

rm(list=ls())

librarian::shelf(tidyverse, here, ggplot2, mvabund)


################################################################################
#set directories and load data
basedir <- "/Volumes/seaotterdb$/kelp_recovery"
figdir <- here::here("analyses","4patch_drivers","Figures")
tabdir <- here::here("analyses","4patch_drivers","Tables")

#load raw dat
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

#load species attribute table
spp_attribute <- readxl::read_excel(file.path(basedir,"data/subtidal_monitoring/raw/spp_attribute_table.xlsx"), sheet = 2)


################################################################################
#process 

#drop species that were never encountered
swath_build1 <- swath_raw %>% dplyr::select(where(~ any(. != 0))) %>%
                        pivot_longer(cols=12:ncol(.), names_to="species", values_to = "counts") %>%
                        mutate(survey_method = "Swath") %>%
                        dplyr::select(survey_method, everything())
upc_build1 <- upc_raw %>% dplyr::select(where(~ any(. != 0))) %>%
                        pivot_longer(cols=12:ncol(.), names_to="species", values_to = "counts") %>%
                        mutate(survey_method = "UPC")%>%
                        dplyr::select(survey_method, everything())
fish_build1 <- fish_raw %>% dplyr::select(where(~ any(. != 0))) %>%
                        pivot_longer(cols=13:ncol(.), names_to="species", values_to = "counts")%>%
                        mutate(survey_method = "Fish")%>%
                        dplyr::select(survey_method, everything())


################################################################################
#find unique species

swath_unique <- swath_build1 %>%
                dplyr::select(survey_method, species) %>% distinct() %>%
                mutate(species = str_to_sentence(gsub("_", " ", species)))

upc_unique <- upc_build1 %>%
  dplyr::select(survey_method, species) %>% distinct() %>%
  mutate(species = str_to_sentence(gsub("_", " ", species)))

fish_unique <- fish_build1 %>%
  dplyr::select(survey_method, species) %>% distinct() %>%
  mutate(species = str_to_sentence(gsub("_", " ", species)))

combined_spp <- rbind(swath_unique, upc_unique, fish_unique)

################################################################################
#process attribute table

spp_tab <- spp_attribute %>% dplyr::select(species = Genusspecies, common_name,
                                           primary_taxonomic, primary_trophic) %>%
                  distinct()%>%
              #drop special characters
              mutate(species = str_to_sentence(gsub("_", " ", species)),
                     species = str_replace_all(species, "-", ""),
                     species = str_replace_all(species, "/", " "),
                     species = str_replace_all(species, "\\(|\\)", ""),
                     common_name = str_to_sentence(gsub("_", " ", common_name)),
                     common_name = str_replace_all(common_name, "-", ""),
                     common_name = ifelse(common_name == "Macrocystis", "Giant kelp",common_name))

              
################################################################################
#join to create full table

spp_table_full <- left_join(combined_spp, spp_tab, by="species")

#identify missing species
spp_table_NA <- spp_table_full %>% filter(is.na(common_name))

#correct names in spp_tab to match what is in database
spp_tab_corrected <- spp_tab %>%
  mutate(species = case_when(
    species == "Cystoseira osmundacea" ~ "Stephanocystis osmundacea",
    species == "Tethya aurantia" ~ "Tethya californiana",
    species == "Pisaster ochraceous" ~ "Pisaster ochraceus",
    #species == "Old_Name3" ~ "Cancridae",
    species == "Lytechinus anamesus" ~ "Lytechinus pictus",
    species == "Crassedoma giganteum" ~ "Crassadoma gigantea",
    species == "Megastrea gibberosa" ~ "Pomaulax gibberosus",
    species == "Loxorhynchus scyra spp" ~ "Loxorhynchus crispatus scyra acutifrons",
    species == "Metridium spp" ~ "Metridium",
    species == "Doriopsilla albopunctata" ~ "Cribrinopsis albopunctata",
    species == "Pugettia spp" ~ "Pugettia foliata",
    species == "Parastichopus parvimensis" ~ "Apostichopus parvimensis",
    species == "Parastichopus californicus" ~ "Apostichopus californicus",
    species == "Balanus nubilis" ~ "Balanus nubilus",
    species == "Pachycerianthus fimbratus" ~ "Pachycerianthus fimbriatus",
    species == "Haliotis wallalensis" ~ "Haliotis walallensis",
    #species == "Aplysia californica" ~ "Aplysia vaccaria",
    species == "Cryptolithoides sitchensis" ~ "Cryptolithodes sitchensis",
    #species == "Old_Name3" ~ "Leptasterias hexactis",
    species == "Urticina spp" ~ "Urticina",
    #species == "Old_Name3" ~ "Phaeophyceae",
    species == "Bryozoan" ~ "Bryozoa",
    species == "Desmarestia spp" ~ "Desmarestia",
    species == "Sponge" ~ "Porifera",
    species == "Phyllospadix spp" ~ "Phyllospadix",
    species == "Green algae" ~ "Chlorophyta",
    species == "Barnacle" ~ "Cirripedia",
    species == "Cucumaria spp" ~ "Cucumaria",
    species == "Tunicate colonial,compund,social" ~ "Tunicate colonial compund social",
    species == "Dictyotales spp" ~ "Dictyoneurum californicum reticulatum",
    #species == "Old_Name3" ~ "Dictyotales",
    species == "Red algae leaflike" ~ "Red algae leaf like",
    species == "Tubeworm mat." ~ "Tubeworm mat",
    species == "Serpulorbis squamigerus" ~ "Thylacodes squamigerus",
    species == "Mussel" ~ "Mytilus spp",
    species == "Sebastes serranoides,flavidus" ~ "Sebastes serranoides flavidus",
    species == "Sebastes chrysomelas" ~ "Sebastes chrysomelas carnatus",
    #species == "Old_Name3" ~ "Sebastes serranoides flavidus melanops",
    species == "Sebastes hopkinsi yoy" ~ "Sebastes hopkinsi",
    species == "Bathymasteridae spp" ~ "Bathymasteridae",
    species == "Sebastes diploproa yoy" ~ "Sebastes diploproa",
    species == "Sebastes dalli" ~ "Sebastes dallii",
    species == "Torpedo californica" ~ "Tetronarce californica",
    species == "Platyrhinoides triseriata" ~ "Platyrhinoidis triseriata",
    species == "Pleuronectidae spp" ~ "Pleuronichthys coenosus",
    species == "Citharichthys spp" ~ "Citharichthys",
    species == "Aplysia californica" ~ "Aplysia californica",
    TRUE ~ species  # Keep the original value if no match
  ))

#rejoin

spp_table_full <- left_join(combined_spp, spp_tab_corrected, by="species") %>%
                    mutate(primary_taxonomic = str_to_sentence(primary_taxonomic),
                           primary_trophic = str_to_sentence(primary_trophic))


write.csv(spp_table_full, file.path(tabdir, "TableS1_spp_table.csv"), row.names = FALSE)













