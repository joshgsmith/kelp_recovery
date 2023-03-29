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

dat <- read_xlsx(file.path(datin, "Forage_data_2016tocurrent.xlsx"))
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










