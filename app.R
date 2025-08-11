#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    https://shiny.posit.co/
#

suppressPackageStartupMessages({
  library(shiny)
  library(SummarizedExperiment)
  library(ggplot2)
  library(dplyr)
  library(Chicago)
  library(cowplot)
  library(stringr)
  library(igraph)
  library(RColorBrewer)
  library(rtracklayer)
  library(netdiffuseR)
  library(IRanges)
  library(tibble)
  library(tidyr)
})

# --- Load Data (CHANGE THESE PATHS) ---
# For deployment, place these files in a 'data' folder and use relative paths
# e.g., file.path("data", "TadBorders.txt")
TadBorders <- read.delim("path/to/your/TadBorders.txt")
Prom <- readRDS("path/to/your/Prom.rds")
Enh <- readRDS("path/to/your/Enh.rds")

# Assuming 'SE', 'E11.5', and 'E10.5' are SummarizedExperiment objects
# You'll need to load these as well.
SE <- readRDS("path/to/your/SE.rds")
E11.5 <- readRDS("path/to/your/E11.5.rds")
E10.5 <- readRDS("path/to/your/E10.5.rds")


# Define UI for application that draws a histogram
ui <- fluidPage(

    # Application title
    titlePanel("Old Faithful Geyser Data"),

    # Sidebar with a slider input for number of bins 
    sidebarLayout(
        sidebarPanel(
            sliderInput("bins",
                        "Number of bins:",
                        min = 1,
                        max = 50,
                        value = 30)
        ),

        # Show a plot of the generated distribution
        mainPanel(
           plotOutput("distPlot")
        )
    )
)

# Define server logic required to draw a histogram
server <- function(input, output) {

    output$distPlot <- renderPlot({
        # generate bins based on input$bins from ui.R
        x    <- faithful[, 2]
        bins <- seq(min(x), max(x), length.out = input$bins + 1)

        # draw the histogram with the specified number of bins
        hist(x, breaks = bins, col = 'darkgray', border = 'white',
             xlab = 'Waiting time to next eruption (in mins)',
             main = 'Histogram of waiting times')
    })
}

# Run the application 
shinyApp(ui = ui, server = server)
