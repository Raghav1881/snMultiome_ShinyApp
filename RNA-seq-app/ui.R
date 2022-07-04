library(shiny)
library(Seurat)
library(ggplot2)
library(uwot)
library(plotly)
library(bslib)

metadata <- colnames(dataset[[]])
feature_list <- c("nCount_RNA", "nCount_SCT", "nFeature_RNA",
                  "nFeature_SCT", "percent.mt", "percent.rpl", "percent.rps")

ui <- fluidPage(theme = bs_theme(version = 4, bootswatch = "minty"),
  titlePanel("scRNA-seq Brain"),
  # Create front main page
  navbarPage(NULL,
    tabPanel(
      HTML("snRNA-seq"),
      sidebarLayout(
        sidebarPanel(
          width = 2,
          selectInput("select", "Cell Categories:", metadata,
          multiple = TRUE),
          checkboxGroupInput("diag", "Diagnoses:",
                              levels(dataset@meta.data[["diagnosis"]])),
          selectInput("gene_input", "Select gene: ", rownames(x = dataset))),
        mainPanel(
          width = 9,
          fluidRow(
            column(
              6, plotOutput("umap_graph")),
            column(
              6, plotOutput("categorical_plot"))),
          fluidRow(
            column(
              6, plotOutput("gene_plot")),
            column(
              6, plotOutput("vln_gene_plot")))))),
    tabPanel("Extra features",
      # Create input for extra features
      fluidRow(
        # Checkbox input for selecting features
        column(
          2, checkboxGroupInput("feats",
                                "Display extra features", feature_list)),
        # Plot features graphs
        column(
          10, plotOutput("features_graph", height = "80vh")))
    )))