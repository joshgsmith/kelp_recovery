

#Joshua G. Smith; jossmith@mbayaq.org

rm(list=ls())


######
#required packages
librarian::shelf(tidyverse, sf, zoo)

#set directories 
basedir <- "/Volumes/seaotterdb$/kelp_recovery/data/census_data"
censusdir <- file.path(basedir,"annual_surveys/processed")
gisdir <- "/Volumes/seaotterdb$/kelp_recovery/data/gis_data"

#read USGS processed data
usgs_orig <- st_read(file.path(censusdir, "usgs_1985_2019_annual_census.shp")) 

#read mpen survey data
mpen_orig <- readRDS(file.path(censusdir, "mpen_annual_counts_2019_23.Rds")) 

#read atos polygons
atos_orig <- st_read(file.path(gisdir, "raw/ATOS/ATOS_polygon_teale83/ATOS_polygon_teale83.shp"))

################################################################################
#workflow

#data preparation
# 1. define area around the mpen surveys
# 2. crop the usgs surveys to the mpen area
# 3. assign ATOS polygons to mpen counts

#moving average
# 4. calculate summary stats at atos level
# 5. prep data for rolling average
# 6. claculate spatial rolling average
# 7. calculate temporal rolling average
# 8. QAQC


####Concerns
# - one mpen atos ID does not match the usgs atos ID. -- the only case is 2019. 
#     mpen data called the location atos 365, but usgs called it 366
# - five mpen points were outside of an ATOS cell - these were joined with the nearest atos cell


################################################################################
# step 1 - define bounding box around mpen survey extent to get the area

#make spatial
mpen_spat <- mpen_orig %>% st_as_sf(coords = c("longitude", "latitude"), crs = 4326)%>%
                            st_transform(crs = 3310)

# Get the original bounding box
bbox_original <- st_bbox(mpen_spat)

# Calculate the width and height of the bounding box
width <- bbox_original$xmax - bbox_original$xmin
height <- bbox_original$ymax - bbox_original$ymin

# Add 1 kilometer padding to the bounding box in CRS 3310
padding <- 1000  
bbox_padded_3310 <- list(xmin = bbox_original$xmin - padding,
                         xmax = bbox_original$xmax + padding,
                         ymin = bbox_original$ymin - padding,
                         ymax = bbox_original$ymax + padding) 

# Convert the bounding box with padding to an sf object in CRS 3310
bbox_padded_3310_sf <- st_bbox(st_sfc(st_polygon(list(rbind(c(bbox_padded_3310$xmin, bbox_padded_3310$ymin),
                                                            c(bbox_padded_3310$xmax, bbox_padded_3310$ymin),
                                                            c(bbox_padded_3310$xmax, bbox_padded_3310$ymax),
                                                            c(bbox_padded_3310$xmin, bbox_padded_3310$ymax),
                                                            c(bbox_padded_3310$xmin, bbox_padded_3310$ymin)))), crs = 3310))


################################################################################
# step 2 - crop the usgs surveys to the mpen area

usgs_mpen <- st_crop(usgs_orig, bbox_padded_3310_sf)

#inspect
plot(usgs_mpen)


################################################################################
# step 3 - assign atos polygons to mpen counts

# note - we have the atos layer, but to be consistent we will use the polygons
# from the usgs historical data

usgs_atos <- usgs_mpen %>% dplyr::select(atos_id, geometry, depth) %>%
              mutate(atos_loc = ifelse(depth == "0 to -30m","onshore",ifelse(
                depth == "-30 to -60m","offshore",atos_id
              )),
              atos_id_loc = paste(atos_id,atos_loc))%>%
              distinct()

#inspect
plot <- ggplot() +
  geom_sf(data = usgs_atos #%>% filter(atos_loc == "offshore")
          ) +
  geom_sf_text(
    data = usgs_atos,
    aes(label = atos_id_loc),
    size = 3,
    nudge_y = 0.02
  ) +
  theme_minimal()
plot


mpen_spat_build1 <- st_join(mpen_spat, usgs_atos)

#inspect
#check if any points are outside a atos polygon
points_outside <- mpen_spat_build1[is.na(mpen_spat_build1$atos_id_loc), ]
#inspect - check for atos ID that do not match between mpen data and usgs
unmatch_atos <- mpen_spat_build1 %>% filter(atos != atos_id)


# Join mpen_spat with usgs_atos using a left join
mpen_spat_build2 <- st_join(mpen_spat, usgs_atos, join = st_nearest_feature)

#clean up
mpen_build3 <- mpen_spat_build2 %>% dplyr::select(-atos,-zone,-hab_id)

################################################################################
# step 4 - calculate the summary stats at atos level


# Group by 'atos_id' and calculate sums for 'ind,' 'smpup,' and 'lgpup'
summarized_data <- mpen_build3 %>%
  st_drop_geometry() %>%
  group_by(year, atos_id_loc) %>%
  summarize(
    sum_ind = sum(ind, na.rm = TRUE),
    sum_smpup = sum(smpup, na.rm = TRUE),
    sum_lgpup = sum(lgpup, na.rm = TRUE)
  )


# Left join summarized_data with usgs_atos_sf by 'atos_id'
mpen_build4 <- left_join(summarized_data, usgs_atos, by = "atos_id_loc")

################################################################################
# step 5 - prepare for rolling average

#clean up and aggregate
mpen_build5 <- mpen_build4 %>%
                mutate(n_indep = sum_ind,
                       n_pup = sum_smpup+sum_lgpup)%>%
                dplyr::select(year, atos_id, atos_loc, n_indep, n_pup)
              


# Get unique years from mpen_build5
unique_years <- unique(mpen_build5$year)

# Replicate rows in usgs_atos for each unique year
usgs_atos_replicated <- usgs_atos %>%
  crossing(year = unique_years)

# join
mpen_build6 <- left_join(usgs_atos_replicated, mpen_build5, by = c("year", "atos_id", "atos_loc")) %>%
                    st_as_sf()

# Replace missing values with 0 for n_indep and n_pup
mpen_build6$n_indep[is.na(mpen_build6$n_indep)] <- 0
mpen_build6$n_pup[is.na(mpen_build6$n_pup)] <- 0


# Inspect
plot(mpen_build6)



################################################################################
# step 6 - calculate rolling average

# Define the window size in terms of atos unit
window_size <- 10 

mpen_build7 <- mpen_build6 %>%
  group_by(year, atos_loc) %>%
  mutate(
    roll_n_indep = zoo::rollapplyr(n_indep, width = 2 * window_size + 1, FUN = mean, fill = NA),
    roll_n_pup = zoo::rollapplyr(n_pup, width = 2 * window_size + 1, FUN = mean, fill = NA)
  ) %>%
  ungroup()


################################################################################
# step 7 - calculate temporal average

#select usgs data to include 
usgs_avg <- usgs_mpen %>% dplyr::filter(year == 2018 | year ==2017) %>%
                #select focal fields
                dplyr::select(year, atos_id, geometry, depth, dens_sm,area) %>%
                mutate(atos_loc = ifelse(depth == "0 to -30m","onshore",ifelse(
                    depth == "-30 to -60m","offshore",atos_id
                  )),
                  n_indep = area*dens_sm,
                  atos_id_loc = paste(atos_id,atos_loc)) %>%
                  dplyr::select(-area, -dens_sm)

#prep mpen for join
mpen_build8 <- mpen_build7 %>% dplyr::select(-n_indep, -n_pup, -roll_n_pup)%>%
                rename(n_indep = roll_n_indep) %>%
                select(year, atos_id, depth, atos_loc, n_indep, 
                        atos_id_loc,geometry
                        )

#join
mpen_build9 <- rbind(mpen_build8, usgs_avg)


# Replace NA with 0 in n_indep
mpen_build9$n_indep[is.na(mpen_build9$n_indep)] <- 0

# Calculate rolling averages
mpen_build10 <- mpen_build9 %>%
  group_by(atos_id_loc) %>%
  arrange(year) %>%
  mutate(
    avg_n_indep = zoo::rollapplyr(n_indep, width = 3, FUN = mean, fill = NA, align = "right")
  ) %>%
  ungroup()




################################################################################
# step 8 - QAQC


mpen_2019 <- mpen_build10 %>% filter(year == 2019)

#select usgs data to include 
usgs_2019 <- usgs_mpen %>% dplyr::filter(year == 2019) %>%
  #select focal fields
  dplyr::select(year, atos_id, geometry, depth, dens_sm,area) %>%
  mutate(atos_loc = ifelse(depth == "0 to -30m","onshore",ifelse(
    depth == "-30 to -60m","offshore",atos_id
  )),
  n_indep = area*dens_sm,
  atos_id_loc = paste(atos_id,atos_loc)) %>%
  dplyr::select(-area, -dens_sm)


View(usgs_2019)
View(mpen_2019)




mpen_total <- mpen_build10 %>% group_by(year) %>% summarize(total = sum(avg_n_indep, na.rm=TRUE))

usgs_total <- usgs_2019 %>% group_by(year) %>% summarize(total = sum(n_indep, na.rm=TRUE))







