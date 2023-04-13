#Joshua G. Smith, jossmith@mbayaq.org

rm(list=ls())

#required packages
librarian::shelf(tidyverse, tidync, sf, rnaturalearth, usethis, raster, tidyterra)

#set directories 
basedir <- "/Volumes/seaotterdb$/kelp_recovery/data"
figdir <- here::here("analyses","figures")

#read landsat dat
landsat_dat <- st_read(file.path(basedir, "kelp_landsat/processed/monterey_peninsula/landsat_mpen_1984_2022_points.shp"))

#read state
ca_counties <- st_read(file.path(basedir, "gis_data/raw/ca_county_boundaries/s7vc7n.shp")) 

#read foraging data
forage_dat <- read.csv(file.path(basedir, "foraging_data/processed/foraging_data_2016_2023.csv")) 

#read forage metadata
forage_meta <- read.csv(file.path(basedir, "foraging_data/processed/forage_metadata.csv")) 


################################################################################
#Step 1 -- process landsat data for join

mpen_kelp <- landsat_dat %>% 
              #filter years that match otter data
              filter(year >= 2016 &
                       #filter data for central coast only
                       latitude >= 36.510140 &
                       latitude <= 36.670574) %>%
  #drop null dat
  filter(!(biomass == "0" | is.na(biomass)))


#transform landsat data to Teale Albers
t_dat <- st_transform(mpen_kelp, crs=3310)

#create grid for raster
r <- terra::rast(t_dat, res=30)  # Builds a blank raster of given dimensions and resolution  

#rasterize
vr <- rasterize(t_dat, r,"biomass", resolution = 30) 

################################################################################
#Step 2 -- process foraging dat for join

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
                              ) %>%
                  filter(!(is.na(lat)))%>%
                  #make spatial
                  st_as_sf(., coords = c("long","lat"), crs=4326)


forage_build4 <- st_transform(forage_build3, crs=3310)


################################################################################
#Step 3 -- calculate distance to nearest point


results <- data.frame(id1 = numeric(), id2 = numeric(), distance = numeric(), year = numeric(), quarter = numeric())


# Loop through each year and quarter
for (year in unique(forage_build4$year)) {
  for (quarter in unique(forage_build4$quarter)) {
    # Get subset of object1 for current year and quarter
    forage_subset <- forage_build4[forage_build4$year == year & forage_build4$quarter == quarter, ]
    # Get subset of object2 for current year and quarter
    kelp_subset <- t_dat[t_dat$year == year & t_dat$quarter == quarter, ]
    # Loop through each point in object1 and find nearest point in object2
    for (i in 1:nrow(forage_subset)) {
     nearest <- st_nearest_feature(forage_subset[i, ], kelp_subset, returnDist=TRUE, units="m")
      
     results <- rbind(results, data.frame( id1 = forage_subset[i, "bout"], id2 = nearest$biomass, 
                                         distance = units::set_units(nearest$dist, m), 
                                          year = year, quarter = quarter))
      
    }
  }
}




################################################################################
#Manual mode

#create empty column 
forage_dist <- forage_build4
forage_dist$nearest_kelp <- NA

dist_fun <- function(y, q) {
  forage_subset <- forage_dist %>% filter(year==y & quarter == q)
  kelp_subset <- t_dat %>% filter(year==y & quarter == q)
  
  #this returns the index (row number) of kelp_subset closest to forage_subset
  out <- st_nearest_feature(forage_subset, kelp_subset)
  #this returns the element-wise distances 
  dist <- st_distance(forage_subset, kelp_subset[out,], by_element=TRUE) %>%
    data.frame() %>% rename("dist" = ".")
  
  df_out <- dist %>% mutate(year = y, q=q,
                                           dist = as.numeric(word(dist,1)),
                                           dist_m = round(dist,2)) %>% dplyr::select(!(dist)) 
  return(df_out)
}


for(yr in unique(forage_dist$year)){
  for(qtr in unique(forage_dist$quarter)){
    dist_out <- dist_fun(yr,qtr)
    forage_dist$nearest_kelp[forage_dist$year==yr & forage_dist$quarter==qtr] <- dist_out$dist_m
  }
}


dist_fun(2019,4)


View(forage_dist)



###trouble shooting
dist_out <- dist_fun(2017,3)

forage_subset <- forage_build4 %>% filter(year==2016 & quarter == 3)
kelp_subset <- t_dat %>% filter(year==2016 & quarter == 3)


dist <- st_distance(forage_subset, kelp_subset[out,], by_element=TRUE, units="m") %>%
            data.frame() %>% rename("dist" = ".")
df_out <- dist %>% mutate(year = 2016, q=3,
                          dist = as.numeric(word(dist,1)),
                          dist_m = round(dist,2)) %>% dplyr::select(!(dist))

forage_dist$nearest_kelp[forage_dist$year==2016 & forage_dist$quarter==3] <- dist_out$dist_m


#this returns the index (row number) of kelp_subset closest to forage_subset
out <- st_nearest_feature(forage_subset, kelp_subset, units="m")
#this returns the element-wise distances 
dist <- st_distance(forage_subset, kelp_subset[out,], by_element=TRUE)

df_out <- as.data.frame(dist) %>% mutate(year = 2016, q=3,
                                         dist = as.numeric(word(dist,1)),
                                         dist_m= round(dist, 2)) %>% select(!(dist))









#######THIS WORKS!!!!!!!!

#create empty column 
forage_dist <- forage_build4


for(yr in unique(forage_dist$year)){
  for(qtr in unique(forage_dist$quarter)){
    #this returns the index (row number) of kelp_subset closest to forage_subset
    out <- st_nearest_feature(forage_dist, t_dat)
    #this returns the element-wise distances 
    dist <- st_distance(forage_dist, t_dat[out,], by_element=TRUE) 
    
    #forage_dist$nearest_kelp[forage_dist$year==yr & forage_dist$quarter==qtr] <- dist$dist_m
    output_df <- cbind(forage_dist, dist)  %>% mutate(
                                                      dist = as.numeric(word(dist,1)),
                                                      dist_m = round(dist,2)) %>% dplyr::select(!(dist)) 
    #join attributes from landsat data
    output_df2 <- cbind(output_df, t_dat[out,])
  }
}






