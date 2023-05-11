#Joshua G. Smith 
#jossmith@mbayaq.org


# Read data
################################################################################

# Clear workspace
rm(list = ls())

librarian::shelf(rerddap, lubridate, tidyverse, beepr, data.table, purrr)


# Directories
basedir <- "/Volumes/seaotterdb$/kelp_recovery"
figdir <- here::here("analyses","4patch_drivers","Figures")
