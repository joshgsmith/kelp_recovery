

#Joshua G. Smith; jossmith@mbayaq.org

rm(list=ls())
#library(pryr)
#pre_memory <- mem_used()  # Capture memory usage before running the script

######
#required packages

library(tidyverse)
library(sf)
library(raster)
library(shiny)
library(tmap)
library(terra)


#set directories 
#setwd(here::here("analyses","2incipient_forests","analyses","foraging_shiny"))

#read landsat dat
#landsat_dat <- st_read("landsat_short.shp") #old
landsat_dat <- st_read("landsat_2016_2023.shp")

#read foraging data
forage_dat <- read.csv("foraging_data_2016_2023.csv")


################################################################################
#Step 1 - build stacked raster

rast_raw <- landsat_dat 

# transform landsat data to Teale Albers
rast_build1 <- st_transform(rast_raw, crs = 3310) %>% 
  mutate(biomass = ifelse(biomass == 0, NA, biomass))  %>% filter (year >= 2016) #reduce years if memory issue

#define blank raster
r <- rast(rast_build1, res=30)

# Get unique quarters and years
quarters <- unique(rast_build1$quarter)
years <- unique(rast_build1$year)

# Loop through each quarter and year
raster_list <- list()
for (y in years) {
  for (q in quarters) {
    # Subset rast_build1 for current quarter and year
    subset <- rast_build1 %>% filter(year == y, quarter == q)
    # Rasterize subset
    rasterized <- rasterize(subset, r, field = "biomass", fun = mean)
    # Add year and quarter as part of the SpatRaster object name
    name <- paste0("year_", y, "_quarter_", q)
    names(rasterized) <- name
    # Add raster to list
    raster_list <- append(raster_list, rasterized)
  }
}

# Stack rasters in list
stacked_raster <- stack(raster_list)

#check names
names(stacked_raster)


################################################################################
#Step 2 - define max kelp extent

#select MPEN as focal region 
plot_dat <- landsat_dat %>% filter ( year == "2019",
                                     latitude >= 36.510140 &
                                       latitude <= 36.670574) 

#define 0 cover as historical kelp footprint
na_dat <- landsat_dat %>% filter (latitude >= 36.510140 &
                                    latitude <= 36.670574, 
                                  biomass == 0)

#transform landsat data to Teale Albers
t_dat <- st_transform(plot_dat, crs=3310) %>% 
  mutate(biomass = ifelse(biomass==0,NA,biomass))%>%
  filter(year == 2019)

kelp_historic <- st_transform(na_dat, crs=3310) # %>% mutate(biomass = ifelse(biomass==0,1,NA))

#create grid
r <- rast(t_dat, res=30)  # Builds a blank raster of given dimensions and resolution  
vr <- rasterize(t_dat, r,"biomass", resolution = 30) 

kelp_na <- rasterize(kelp_historic, r,"biomass", resolution = 30) 

################################################################################
#Step 3 -- process foraging data

forage_build1 <- forage_dat %>%
  filter (year >= 2016) %>%
  mutate(quarter = ifelse(month == 1 | month == 2 | month == 3, 1,
                          ifelse(month == 4 | month == 5 | month == 6,2,
                                 ifelse(month == 7 | month == 8 | month == 9,3,4)))) # %>%
  #filter urchin prey only
 # filter(prey == "pur" | prey == "mus") 

focal_patch <-  forage_dat %>%
  filter (year >=2016)%>%
  mutate(quarter = ifelse(month == 1:3,1,
                          ifelse(month == 4:6,2,
                                 ifelse(month==7:9,3,4)))) %>%
  #filter urchin prey only
 # filter(prey == "pur" | prey == "mus") %>%
  #define focal patch
  group_by(bout) %>%
  summarize(n_dives = n()) %>%
  mutate(focal_patch = ifelse(n_dives > 3, "yes","no"))


#join

forage_build2 <- left_join(forage_build1, focal_patch, by="bout")


#aggregate data by bout

forage_build3 <- forage_build2 %>%
  group_by(year, quarter, bout, focal_patch, prey) %>%
  summarize(n_prey = sum(preynum),
            lat = mean(lat),
            long = mean(long)
  ) %>%
  filter(!(is.na(lat)))%>%
  #make spatial
  st_as_sf(.,coords = c("long","lat"), crs=4326)


forage_plot_dat <- forage_build3

################################################################################
#Build shiny

ui <- fluidPage(
  tags$head(
    tags$style(
      HTML(".leaflet-container { height: 600px !important; }"),
      HTML("HTML, body {background-color: transparent;}")
    )
  ),
  fluidRow(
    column(width = 4,
           sliderInput("year_toggle", "Filter by year:", min = min(forage_plot_dat$year), max = max(forage_plot_dat$year), value = min(forage_plot_dat$year), step = 1, sep = "", width = "90%"),
           sliderInput("quarter_toggle", "Filter by quarter:", min = 1, max = 4, value = 1, step = 1, sep = "", width = "90%"),
           selectInput("prey_toggle", "Filter by prey:", c("All", unique(forage_plot_dat$prey)), multiple = TRUE, selected = "All"),
           selectInput("focal_patch_toggle", "Filter by focal patch:", c("Yes", "No", "N/A"), selected = "N/A")
    ),
    column(width = 8,
           tmap::tmapOutput("map"),
           br(),
           p(tags$b("FOR INTERNAL USE ONLY")),
           p("Created by: Joshua G. Smith - Monterey Bay Aquarium"),
           p("Source files:"),
           p("Santa Barbara Coastal LTER, T. Bell, K. Cavanaugh, and D. Siegel. 2023. SBC LTER: Time series of quarterly NetCDF files of kelp biomass in the canopy from Landsat 5, 7 and 8, since 1984 (ongoing) ver 20. Environmental Data Initiative. https://doi.org/10.6073/pasta/41f330ccf66fa8c05fc851862e69b1da (Accessed 2023-06-02)."),
           p("Sea otter foraging - Monterey Bay Aquarium"),
    )
  )
)



server <- function(input, output, session) {
  
  filtered_data <- reactive({
    if (input$focal_patch_toggle == "Yes") {
      filter(forage_plot_dat, focal_patch == "yes")
    } else if (input$focal_patch_toggle == "No") {
      filter(forage_plot_dat, focal_patch == "no")
    } else {
      forage_plot_dat
    }
  })
  
  filtered_data_year_quarter <- reactive({
    filter(filtered_data(), year == input$year_toggle, quarter == input$quarter_toggle)
  })
  
  filtered_data_prey <- reactive({
    if ("All" %in% input$prey_toggle) {
      filtered_data_year_quarter()
    } else {
      filter(filtered_data_year_quarter(), prey %in% input$prey_toggle)
    }
  })
  
  tm_map_reactive <- eventReactive(list(input$focal_patch_toggle, input$year_toggle, input$quarter_toggle, input$prey_toggle), {
    
    # Check if the selected year and quarter combination is available in stacked_raster
    stacked_raster_layer_name <- paste0("year_", input$year_toggle, "_quarter_", input$quarter_toggle)
    if (!stacked_raster_layer_name %in% names(stacked_raster)) {
      tm_map <- tm_shape(kelp_na) + 
        tm_raster(palette = "Blues", title = "Kelp footprint \n(all time max)", labels = "")
    } else if (nrow(filtered_data_prey()) > 0) {
      tm_map <- tm_shape(kelp_na) +
        tm_raster(palette = "Blues", title = "Kelp footprint \n(all time max)", labels = "") +
        tm_shape(stacked_raster[[stacked_raster_layer_name]]) + 
        tm_raster(breaks = c(0, 1000, 2000, 4000, 8000), title = "Kelp biomass (Kg)") +
        tm_shape(filtered_data_prey()) +
        tm_symbols(scale = 0.1, col = "black")
    } else {
      tm_map <- tm_shape(kelp_na) +
        tm_raster(palette = "Blues", title = "Kelp footprint \n(all time max)", labels = "") +
        tm_shape(stacked_raster[[stacked_raster_layer_name]]) + 
        tm_raster(breaks = c(0, 1000, 2000, 4000, 8000), title = "Kelp biomass (Kg)")
    }
    
    tm_map
  })
  
  
  output$map <- renderTmap({
    tm_map_reactive()
  })
}

shinyApp(ui, server)

 
#post_memory <- mem_used()  # Capture memory usage after running the script
#memory_used <- post_memory - pre_memory
