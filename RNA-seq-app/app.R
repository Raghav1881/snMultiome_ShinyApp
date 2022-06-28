library(shiny)
library(Seurat)
library(ggplot2)
library(uwot)
library(plotly)

metadata <- colnames(pbmc_tutorial[[]])

ui <- fluidPage(
  titlePanel("scRNA-seq Mouse Brain"),
  navbarPage(NULL,
    tabPanel("Sample selector",
      column(3,
        fluidRow(selectInput("select", "Catagorical Variables:", metadata, 
                             multiple = TRUE)),
        fluidRow(checkboxGroupInput("diag", "Diagnoses:", levels(x = pbmc_tutorial)))
      ),
      column(6,
      )
    )
  ))

server <- function(input, output) {

}

shinyApp(ui = ui, server = server)
