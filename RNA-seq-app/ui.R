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
    tabPanel("Sample selector",
      # Create inputs for selecting categorical variables and the diagnoses
      column(12,
        fluidRow(selectInput("select", "Categorical Variables:", metadata,
                             multiple = TRUE)),
        fluidRow(checkboxGroupInput("diag", "Diagnoses:",
                                    levels(dataset@meta.data[["diagnosis"]])))
      ),
      fluidRow(
        # Create input for extra features
        column(2,
          checkboxGroupInput("feats", "Display extra features",
                            feature_list)),
        # Plot features graphs
        column(10,
          plotOutput("features_graph"))
      )
    )
  ))
