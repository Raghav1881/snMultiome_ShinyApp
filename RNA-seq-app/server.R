library(shiny)
library(Seurat)
library(ggplot2)
library(uwot)
library(plotly)
dataset <- readRDS("30k_snRNA.RDS")

server <- function(input, output) {
    if ((length(observe(input$feats))) > 1) {
      output$features_graph <- renderPlot({
        VlnPlot(dataset, features = input$feats,
                pt.size = 0, ncol = length(input$feats)) +
            theme(axis.title.x = element_blank())
      })
    }
}