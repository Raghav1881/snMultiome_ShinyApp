library(shiny)
library(Seurat)
library(ggplot2)
library(uwot)
library(plotly)

gene_list <- rownames(x = pbmc_tutorial)

ui <- fluidPage(
  titlePanel("Seurat Gene Selectornator"),
  navbarPage(NULL,
    tabPanel("UMAP Plot", 
      column(3,
        selectInput("geneid", "Gene:", choices = gene_list)
      ), 
      column(9,
        plotOutput("umap_graph"))
      ),
    tabPanel("Violin Plot", 
      column(3,
        selectInput("geneid", "Gene:", choices = gene_list)
      ), 
      column(9,
        plotOutput("vln_plot"))
      ),
    tabPanel("Feature Plot",
      column(3,
        selectInput("geneid", "Gene:", choices = gene_list)
      ), 
      column(9,
        plotOutput("ftr_plot"))),
    tabPanel("Coexpression Plot", 
      fluidRow(
        plotOutput("coexpress_plot")
      ),
      fluidRow(
        selectInput("gene1", "Gene 1:", choices = gene_list), 
        selectInput("gene2", "Gene 2:", choices = gene_list) 
      )
    )
  ))

server <- function(input, output) {
  output$umap_graph <- renderPlot({
    plot <- FeaturePlot(pbmc_tutorial, features = input$geneid)
    HoverLocator(plot = plot, 
                 information = FetchData(pbmc_tutorial, 
                                         vars = c("ident", "PC_1", "nFeature_RNA")))
  })
  
  output$vln_plot <- renderPlot({
    VlnPlot(pbmc_tutorial, features = input$geneid, pt.size = 0)
  })
  
  output$ftr_plot <- renderPlot({
    plot <- FeaturePlot(pbmc_tutorial, features = input$geneid)
    LabelClusters(plot = plot, id = "ident")
  })
  
  output$coexpress_plot <- renderPlot({
    FeaturePlot(pbmc_tutorial, features = c(input$gene1, input$gene2), 
                blend = TRUE)
  })
}

# Run the application 
shinyApp(ui = ui, server = server)
