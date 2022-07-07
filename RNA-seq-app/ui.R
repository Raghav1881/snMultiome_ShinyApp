library(shiny)
library(Seurat)
library(ggplot2)
library(uwot)
library(plotly)
library(bslib)

diagnosis_list <- dataset@meta.data[["diagnosis"]]
metadata <- colnames(dataset[[]])
feature_list <- c("nCount_RNA", "nCount_SCT", "nFeature_RNA",
                  "nFeature_SCT", "percent.mt", "percent.rpl", "percent.rps")
gene_list <- rownames(x = dataset)

ui <- fluidPage(theme = bs_theme(version = 4, bootswatch = "minty"),
  titlePanel("scRNA-seq Brain"),
  # Create front main page
  navbarPage(
    NULL, theme = bs_theme(version = 4, bootswatch = "minty"),
    tabPanel("snRNA-seq",
      HTML("snRNA-seq UMAP Data"),
      sidebarLayout(
        sidebarPanel(
          width = 3,
          selectizeInput("select", "Cell Categories:",
                        choices = NULL,
                        multiple = TRUE),
          checkboxGroupInput("diag", "Diagnoses:", levels(diagnosis_list)),
          selectizeInput("gene_input", "Select gene:", choices = NULL)),
        mainPanel(
          width = 9,
          fluidRow(
            column(
              6, plotOutput("umap_graph")),
            column(
              6, plotOutput("categorical_plot"))),
          fluidRow(
            column(
              6, plotOutput("gene_plot")))))),
    # Seperate tab for comparing diag vs gene expr data
    tabPanel("Diagnosis & Gene Expr",
      HTML("Diagnosis and Gene Expression Analysis"),
      sidebarLayout(
        sidebarPanel(
          width = 3,
          selectizeInput("genediag", "Select gene", choices = NULL),
          checkboxGroupInput("diagchk", "Select diagnosis",
                      levels(diagnosis_list)),
          checkboxGroupInput("celltype", "Select cell types",
                            levels(dataset))),
        # Generate violin plot output based on gene, diagnosis, and cell types
        mainPanel(
          width = 9, plotOutput("vln_gene_plot")))),
    tabPanel("Extra features",
      # Create input for extra features
      fluidRow(
        # Checkbox input for selecting features
        column(
          2, checkboxGroupInput("feats",
                                "Display extra features", feature_list)),
        # Plot features graphs
        column(
          10, plotOutput("features_graph",
          height = "80vh",
          width = "80vh")))
    )))