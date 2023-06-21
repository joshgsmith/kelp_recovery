
#Joshua G. Smith; jossmith@mbayaq.org

rm(list=ls())

######
#required packages

library(tidyverse)
library(sf)
library(raster)
library(shiny)
library(tmap)
library(terra)

#set directories 
setwd(here::here("analyses","2forage_behavior","analyses","foraging_shiny"))

#read landsat dat
#landsat_dat <- st_read("landsat_short.shp") #old
landsat_dat <- st_read("landsat_2016_2023.shp")

ca_counties <- st_read("/Volumes/seaotterdb$/kelp_recovery/data/gis_data/raw/ca_county_boundaries/s7vc7n.shp") 


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
library(shiny)
library(leaflet)
library(raster)

# Load the stacked raster


# Extract year and quarter information from layer names
layer_names <- names(stacked_raster)
layer_info <- strsplit(layer_names, "_")
years <- unique(sapply(layer_info, function(x) as.integer(x[2])))
quarters <- unique(sapply(layer_info, function(x) as.integer(x[4])))

# Define the UI
ui <- fluidPage(
  titlePanel("Stacked Raster Explorer"),
  sidebarLayout(
    sidebarPanel(
      sliderInput("year", "Year", min(years), max(years), min(years), step = 1),
      sliderInput("quarter", "Quarter", min(quarters), max(quarters), min(quarters), step = 1)
    ),
    mainPanel(
      leafletOutput("map")
    )
  )
)

# Define the server
server <- function(input, output) {
  
  # Reactive expression to filter the stacked raster based on year and quarter
  selected_layer <- reactive({
    layer_name <- paste0("year_", input$year, "_quarter_", input$quarter)
    layer_index <- grep(layer_name, layer_names)
    return(stacked_raster[[layer_index]])
  })
  
  # Render the leaflet map
  output$map <- renderLeaflet({
    leaflet() %>%
      addProviderTiles("Esri.WorldImagery") %>%
      addRasterImage(selected_layer(), opacity = 0.7)
  })
}

# Run the app
shinyApp(ui, server)




################################################################################
#test with GPT === PlotLy thingy
if (!require(plotly)) {
  install.packages("plotly")
}

# Extract year and quarter information from layer names
layer_names <- names(stacked_raster)
layer_info <- strsplit(layer_names, "_")
years <- unique(sapply(layer_info, function(x) as.integer(x[2])))
quarters <- unique(sapply(layer_info, function(x) as.integer(x[4])))

# Define the UI
ui <- fluidPage(
  titlePanel("Stacked Raster Explorer"),
  sidebarLayout(
    sidebarPanel(
      sliderInput("year", "Year", min(years), max(years), min(years), step = 1),
      sliderInput("quarter", "Quarter", min(quarters), max(quarters), min(quarters), step = 1)
    ),
    mainPanel(
      plotlyOutput("plot", height = "600px")
    )
  )
)

# Define the server
server <- function(input, output) {
  
  # Reactive expression to filter the stacked raster based on year and quarter
  selected_layer <- reactive({
    layer_name <- paste0("year_", input$year, "_quarter_", input$quarter)
    layer_index <- grep(layer_name, layer_names)
    return(stacked_raster[[layer_index]])
  })
  
  # Render the 3D plot
  output$plot <- renderPlotly({
    heightmap <- raster::as.matrix(selected_layer())
    plot_map(heightmap)
  })
}

# Function to plot the 3D map using plotly
plot_map <- function(heightmap) {
  library(plotly)
  
  x <- 1:ncol(heightmap)
  y <- 1:nrow(heightmap)
  z <- heightmap
  
  # Create 3D surface plot
  plot_ly(x = x, y = y, z = z, type = "surface") %>%
    layout(scene = list(
      camera = list(
        eye = list(x = -2, y = -2, z = 1)
      )
    ))
}

# Run the app
shinyApp(ui, server)


################################################################################
#test with GPT === .KML


# Install required packages if not already installed
# install.packages(c("raster", "rasterVis", "XML"))

library(raster)
library(rasterVis)
library(XML)

# Load the stacked raster file
#stacked_raster 

# Extract layer names and corresponding metadata
layer_names <- names(stacked_raster)
metadata <- strsplit(layer_names, "_")

# Create a new KML document
kml <- newXMLDoc()

# Create a KML root node
kml_root <- newXMLNode("kml", parent = kml)

# Add the XML namespace
XML::addNamespace(kml_root, "http://www.opengis.net/kml/2.2")

# Create a document node
doc_node <- newXMLNode("Document", parent = kml_root)

# Loop through each layer and create a KML overlay for it
for (i in 1:length(layer_names)) {
  layer <- stacked_raster[[i]]
  layer_name <- layer_names[i]
  layer_metadata <- metadata[[i]]
  
  # Extract year and quarter from metadata
  year <- layer_metadata[grep("year", layer_metadata)][[1]]
  quarter <- layer_metadata[grep("quarter", layer_metadata)][[1]]
  
  # Create a KML GroundOverlay for the layer
  overlay_node <- newXMLNode("GroundOverlay", parent = doc_node)
  
  # Create a name node
  name_node <- newXMLNode("name", parent = overlay_node)
  xmlValue(name_node) <- layer_name
  
  # Create an icon node
  icon_node <- newXMLNode("Icon", parent = overlay_node)
  
  # Create a href node (path to the raster layer)
  href_node <- newXMLNode("href", parent = icon_node)
  xmlValue(href_node) <- paste0("path/to/", layer_name, ".tif")
  
  # Create a LatLonBox node
  lat_lon_box_node <- newXMLNode("LatLonBox", parent = overlay_node)
  
  # Add coordinates of the raster extent
  extent <- extent(layer)
  south <- extent[["ymin"]]
  north <- extent[["ymax"]]
  west <- extent[["xmin"]]
  east <- extent[["xmax"]]
  
  newXMLNode("north", parent = lat_lon_box_node, value = north)
  newXMLNode("south", parent = lat_lon_box_node, value = south)
  newXMLNode("east", parent = lat_lon_box_node, value = east)
  newXMLNode("west", parent = lat_lon_box_node, value = west)
  
  # Create a description node with the metadata
  description <- paste("Year:", year, "Quarter:", quarter)
  description_node <- newXMLNode("description", parent = overlay_node)
  xmlValue(description_node) <- description
}

# Save the KML file
saveXML(kml, file = "path/to/output.kml")



################################################################################
#test with gganimate

library(ggplot2)
library(gganimate)
library(sf)


# Extract the metadata from layer names
layer_names <- names(stacked_raster)
metadata <- strsplit(layer_names, "_")
metadata <- do.call(rbind, metadata)
metadata <- as.data.frame(metadata)
colnames(metadata) <- c("prefix", "year", "quarter")

# Check if the number of layers and metadata rows match
if (length(layer_names) != nrow(metadata)) {
  stop("Number of layers in the raster stack and metadata rows do not match.")
}

# Convert the stacked raster to a data frame with coordinates
raster_df <- raster::as.data.frame(stacked_raster, xy = TRUE)

# Create a layer column in the raster data frame
raster_df$layer <- layer_names[rep(1:length(layer_names), each = nrow(raster_df) / length(layer_names))]

# Merge the metadata with the raster data frame
merged_df <- merge(raster_df, metadata)

# Convert the data frame to an sf object
raster_sf <- st_as_sf(merged_df, coords = c("x", "y"))

# Create the initial plot with ggplot
p <- ggplot() +
  geom_sf(data = raster_sf, aes(fill = layer))

# Add animation settings and render the animation
animation <- p +
  transition_states(
    year,
    transition_length = 2,
    state_length = 1
  ) +
  enter_fade() +
  exit_fade()

animate(animation)
