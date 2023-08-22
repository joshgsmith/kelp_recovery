#Joshua G. Smith, jossmith@mbayaq.org

rm(list=ls())


# install.packages("librarian")
librarian::shelf(tidyverse, tidync, sf, rnaturalearth, usethis, raster)

#set directories 
datin <- "/Volumes/seaotterdb$/kelp_recovery/data/kelp_landsat/processed/monterey_peninsula"
datout <- here::here("analyses","2forage_behavior","analyses","foraging_shiny")
figdir <- here::here("analyses","figures")

landsat_dat <- st_read(file.path(datin, "landsat_mpen_1984_2023_points_withNAs.shp"))


################################################################################
#process landsat for focal regions and years

landsat_2016_2023 <- landsat_dat %>% filter(year >= 2016,
                                            latitude >= 36.510140 &
                                              latitude <= 36.670574) 


st_write(landsat_2016_2023, file.path(datout,"landsat_2016_2023.shp"))
