library(Seurat)
library(shinyWidgets)
library(shiny)
library(bslib)

mtheme <- bs_theme(version = 4, bootswatch = "minty")
feature_list <- c("nCount_RNA", "nCount_SCT", "nFeature_RNA",
                  "nFeature_SCT", "percent.mt", "percent.rpl", "percent.rps")

ui <- fluidPage(theme = mtheme,
  titlePanel(
    "snATAC- and snRNA-Seq Atlas of ALS/FTLD Orbitofrontal Cortex"),
    navbarPage(
      NULL, theme = mtheme,
      tabPanel("UMAP-RNA"),
      tabPanel("UMAP-ATAC"),
      tabPanel("Coverage/Violin Plots",
        HTML("Diagnosis and Gene Expression Analysis"),
        sidebarLayout(
          sidebarPanel(
            width = 3,
            selectizeInput("genediag",
                        "Select gene",
                        choices = NULL),
            selectizeInput("diagchk",
                        "Select diagnosis",
                        choices = NULL),
            checkboxGroupInput("celltype", "Select cell types",
                                levels(dataset$celltype),
                                selected = "Oligodendrocytes")),
          mainPanel(
            width = 9,
            fluidRow(
              verbatimTextOutput("test"), verbatimTextOutput("test2")),
            fluidRow(
              column(
                6, plotOutput("violin1")),
              column(
                6, plotOutput("violin2")))))),
      tabPanel("Quality Control",
        HTML("Cell Quality Control Metrics"),
        # Create input for extra features
        fluidRow(
          sidebarPanel(
            # Checkbox input for selecting features
            width = 3, checkboxGroupInput("feats", "Display extra features",
                                          feature_list)),
          # Plot features graphs
          mainPanel(
            width = 9, plotOutput("features_graph",
            height = "80vh"))))
    )
)