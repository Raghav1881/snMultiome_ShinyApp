library(shiny)
library(ggplot2)
library(Seurat)
library(uwot)
library(shinythemes)

files.dir <- "~/Documents/GitHub/LMP-RShiny-App/output/"

#Open file location and create a dropdown menu with file options in output dir
rdsfile <- list.files(files.dir, pattern = "\\.rds$")

ui <- fluidPage(theme = shinytheme("sandstone"), pageWithSidebar(
  headerPanel("Seurat Objs"),
  sidebarPanel(
  #Choose dropdown items
    selectInput("dataset", "Choose a dataset:",
                choices = rdsfile),
  ),
  mainPanel(
    tabsetPanel(
      #Individual tabs for each graph type
      tabPanel("Cell types", plotOutput("gene")),
      tabPanel("Individuals", plotOutput("feats"))
    )
  )
))

server <- shinyServer(function(input, output) {
  #Create reactive obj with Seurat obj
  datasetInput <- reactive({
    df <- readRDS(paste0(files.dir, input$dataset))
    return(df)
  })
  #Generate a plot for features plot
  output$feats <- renderPlot({
    dataset <- datasetInput()
    FeaturePlot(dataset, features = c("MS4A1", "GNLY", "CD3E", "CD14", "FCER1A", "FCGR3A", "LYZ", "PPBP",
                                      "CD8A"))
  })
  #Generate a plot for all clusters with identities
  output$gene <- renderPlot({
    dataset <- datasetInput()
    DimPlot(dataset, reduction = "umap")
  })
})

shinyApp(ui = ui, server = server)