#Joshua G. Smith; jossmith@mbayaq.org

rm(list=ls())

######
######
#required packages
librarian::shelf(tidyverse, terra, sf, units, smoothr)

#set directories 
basedir <- "/Volumes/seaotterdb$/kelp_recovery/data"
output <- here::here("analyses","2incipient_forests","output")

#see
#https://cran.r-project.org/web/packages/smoothr/vignettes/smoothr.html