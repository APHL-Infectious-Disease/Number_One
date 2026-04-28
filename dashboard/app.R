#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    https://shiny.posit.co/
#

library(shiny)
library(htmltools)
library(tidyverse)

# Define UI for application that draws a barplot
ui <- fluidPage(
  
  # Application title
  titlePanel("Enteric Virus Reads"),
  
  #UI input - sliders (date), button, arranging. - V 
  
  # Sidebar with a slider input for number of bins 
  sidebarLayout(
    sidebarPanel(
      selectInput("fill_var",
                  "Color By:",
                  choices = c("scientific_name", "sample_srx", "year"))
    ),
    
    # Show a plot of the generated distribution
    mainPanel(
      plotOutput("barPlot")
    )
  )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
  
  # Load data
  kraken_summary <- read.csv("results/summary/kraken_summary-2026-04-28.csv")
  
  output$barPlot <- renderPlot({
    
    # Parse input values (to deal with non-standard ggplot aes nomenclature)
    #color_by <- input$color_by
    
    # generate bins based on input$bins from ui.R
    #bins <- seq(min(x), max(x), length.out = input$bins + 1)
    
    # draw the histogram with the specified number of bins
    kraken_summary %>%
      # Drop unclassified from results
      mutate(year = year(as.Date(summary_date))) %>% #using summary date just to test, but this needs to be changed to the collection date from the metadata
      filter(scientific_name != "unclassified") %>%
      ggplot() +
      geom_bar(aes_string(x = "sample_srx", y="count", fill = input$fill_var), stat = "identity") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      labs(
        y = "read count"
      )
  })
  
  #import data and transform - Script from Hanley
  # tabulated data - data table on a tab
  # ggplot - taxon breakdown by state AND National, and frequency
  # ggplot - time bound
  #
  
  
  
  
  # mapping - Leaflet map. V 
  # state shading? organism
  
  
  
}

#save html
save_html(ui, file = "enteric_ww_dash.html")
# Run the application 
shinyApp(ui = ui, server = server)
