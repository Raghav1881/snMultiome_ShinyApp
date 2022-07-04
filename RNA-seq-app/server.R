library(shiny)
library(Seurat)
library(ggplot2)
library(uwot)
library(plotly)

server <- function(input, output) {

  # Generate output for features plots
  output$features_graph <- renderPlot({
    VlnPlot(dataset, features = input$feats,
            pt.size = 0, ncol = 3) &
        theme(axis.title.x = element_blank())
  })

  output$umap_graph <- renderPlot({
    DimPlot(dataset)
  })

  output$categorical_plot <- renderPlot({
    UMAPPlot(dataset, group.by = input$select)
  })

  output$gene_plot <- renderPlot({
    FeaturePlot(dataset, features = input$gene_input)
  })

  output$vln_gene_plot <- renderPlot({
    VlnPlot(dataset, features = input$gene_input,
            idents = c(input$diag)) &
        theme(axis.title.x = element_blank())
  })
}