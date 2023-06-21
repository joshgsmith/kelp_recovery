library(shiny)

# Define the UI
ui <- fluidPage(
  titlePanel("Satellite Imagery Explorer"),
  mainPanel(
    tags$iframe(src = "map.html", style = "height: 500px; width: 100%; border: none;")
  )
)

# Run the app
shinyApp(ui, function(input, output, session) {})