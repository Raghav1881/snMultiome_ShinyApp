library(shiny)
library(Seurat)
library(ggplot2)
library(uwot)
library(plotly)

metadata <- colnames(dataset[[]])
feature_list <- c("nCount_RNA", "nCount_SCT", "nFeature_RNA",
                  "nFeature_SCT", "percent.mt", "percent.rpl", "percent.rps")

ui <- fluidPage(
  titlePanel("scRNA-seq Mouse Brain"),
  # Create front main page
  navbarPage(NULL,
    tabPanel("snRNA-seq",
      sidebarLayout(
        sidebarPanel(
          selectInput("select", "Categorical Variables:", metadata,
          multiple = TRUE),
          checkboxGroupInput("diag", "Diagnoses:",
                              levels(dataset@meta.data[["diagnosis"]])),
          selectInput("gene_input", "Select gene: ", rownames(x = dataset))),
        mainPanel(
          fluidRow(plotOutput("umap_graph"))))),
    tabPanel("Extra features",
      # Create input for extra features
      column(2,
        checkboxGroupInput("feats", "Display extra features",
                          feature_list)),
      # Plot features graphs
      column(10,
        plotOutput("features_graph")))
    ))