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
      tabPanel("UMAP-RNA",
        HTML("Dimensional and Feature Plot of snRNA-seq Data"),
        sidebarLayout(
          sidebarPanel(
            width = 3,
            selectizeInput("genediag1",
                        "Select gene",
                        choices = NULL),
            selectizeInput("diagchk1",
                        "Select diagnosis",
                        choices = NULL),
            checkboxGroupInput("celltype1", "Select cell types",
                                levels(dataset$celltype),
                                selected = "Oligodendrocytes")),
          mainPanel(
            width = 9,
            fluidRow(
              column(
                6,
                plotOutput("dimPlotRNA")),
              column(
                6,
                plotOutput("featPlotRNA"))),
            fluidRow(
              column(
                6,
                plotOutput("dimPlotRNACtrl")),
              column(
                6,
                plotOutput("dimPlotRNADiag")))
        ))),
      tabPanel("UMAP-ATAC"),
      tabPanel("Coverage/Violin Plots",
        HTML("Diagnosis and Gene Expression Analysis"),
        sidebarLayout(
          sidebarPanel(
            width = 3,
            selectizeInput("genediag2", "Select gene",
                            choices = NULL),
            selectizeInput("diagchk2","Select diagnosis",
                            choices = NULL),
            checkboxGroupInput("celltype2", "Select cell types",
                                levels(dataset$celltype),
                                selected = "Oligodendrocytes")),
          mainPanel(
            width = 9,
            fluidRow(
              verbatimTextOutput("test"), verbatimTextOutput("test2")),
            fluidRow(
              column(
                6, plotOutput("violin1", height = "40vh")),
              column(
                6, plotOutput("violin2", height = "40vh")))))),
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
            height = "60vh"))))
    )
)