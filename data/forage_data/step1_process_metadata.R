#Joshua G. Smith
#jossmith@mbayaq.org
#Script initiated March 29, 2023

rm(list=ls())

#load packages
librarian::shelf(tidyverse, readxl, here, janitor)

#set directories 
datain <- "/Volumes/seaotterdb$/kelp_recovery/data/foraging_data/raw"
dataout <- "/Volumes/seaotterdb$/kelp_recovery/data/foraging_data/processed"

#get metadata 
meta_raw <- read_xlsx(file.path(datin, "ForageKey_lookup_table.xlsx"))


################################################################################
#step 1 - tidy data


meta_build1 <- meta_raw %>%
                janitor::clean_names() %>%
                #Sentence case
                mutate(category_code = ifelse(data_type == "Prey_vars",
                                          str_to_sentence(category_code), category_code),
                       category_long = str_to_sentence(category_long),
                       description = str_to_sentence(description))

################################################################################
#write clean metadata to .csv

write.csv(meta_build1, file.path(dataout, "forage_metadata.csv"),row.names=F)
                






