library(shiny)
library(Seurat)
library(ggplot2)
library(uwot)

gene_list <- rownames(x = pbmc_tutorial)
# Define UI for application that draws a histogram
ui <- fluidPage(
  headerPanel("Seurat Gene Selectornator"),
  sidebarPanel(
    selectInput("geneid", "Gene:", choices = gene_list)
  ),
  mainPanel(
    tabsetPanel(
      tabPanel("UMAP Plot", plotOutput("umap_graph")),
      tabPanel("Violin Plot", plotOutput("vln_plot")),
      tabPanel("Feature Plot", plotOutput("ftr_plot"))
    ))
    
)

# Define server logic required to draw a histogram
server <- function(input, output) {
  output$umap_graph <- renderPlot({
    UMAPPlot(pbmc_tutorial)
  })
  
  output$vln_plot <- renderPlot({
    VlnPlot(pbmc_tutorial, features = input$geneid)
  })
  
  output$ftr_plot <- renderPlot({
    FeaturePlot(pbmc_tutorial, features = input$geneid)
  })
}

# Run the application 
shinyApp(ui = ui, server = server)
