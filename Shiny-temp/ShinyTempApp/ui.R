library(Seurat)
library(shinyWidgets)
library(shiny)
library(bslib)
library(shinycssloaders)
library(shinycustomloader)

mtheme <- bs_theme(version = 5, bootswatch = "flatly")
feature_list <- c("nCount_RNA", "nCount_SCT", "nFeature_RNA",
                  "nFeature_SCT", "percent.mt", "percent.rpl", "percent.rps")

ui <- fluidPage(theme = mtheme,
  titlePanel(
    "snATAC- and snRNA-Seq Atlas of ALS/FTLD Orbitofrontal Cortex"),
    navbarPage(
      NULL,
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
            checkboxGroupInput("celltype1",
                               "Select cell types",
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
                12,
                plotOutput("dimPlotRNACtrl")))
        ))),
      tabPanel("UMAP-ATAC",
        HTML("Dimensional and Feature Plot of snATAC-seq Data"),
        sidebarLayout(
          sidebarPanel(
            width = 3,
            selectizeInput("genediag2",
                           "Select gene",
                           choices = NULL),
            selectizeInput("diagchk2",
                           "Select diagnosis",
                           choices = NULL),
            checkboxGroupInput("celltype2", "Select cell types",
                                levels(dataset$celltype),
                                selected = "Oligodendrocytes")),
          mainPanel(
            width = 9,
            fluidRow(
              column(
                6,
                plotOutput("dimPlotATAC")),
              column(
                6,
                plotOutput("featPlotATAC"))),
            fluidRow(
              column(
                12,
                plotOutput("dimPlotATACCtrl")))

        ))),
      tabPanel("Coverage/Violin Plots",
        HTML("Diagnosis and Gene Expression Analysis"),
        sidebarLayout(
          sidebarPanel(
            width = 3,
            selectizeInput("genediag3", "Select gene",
                            choices = NULL),
            selectizeInput("diagchk3","Select diagnosis",
                            choices = NULL),
            checkboxGroupInput("celltype3", "Select cell types",
                                levels(dataset$celltype),
                                selected = "Oligodendrocytes")),
          mainPanel(
            width = 9,
            h1("Coverage Plot"),
            fluidRow(
              plotOutput("coverage_plot",
                          height = "60vh"),
              br(),
              hr()),
            fluidRow(
              h1("Violin Plots"),
              br(),
              column(
                6, plotOutput("violin1",
                              height = "40vh")),
              column(
                6, plotOutput("violin2",
                              height = "40vh")))
            )
        )),
      tabPanel("Quality Control",
        HTML("Cell Quality Control Metrics"),
        # Create input for extra features
        fluidRow(
          sidebarPanel(
            # Checkbox input for selecting features
            width = 3, checkboxGroupInput("feats",
                                          "Display extra features",
                                          feature_list)),
          # Plot features graphs
          mainPanel(
            width = 9,
            plotOutput("features_graph",
                        height = "60vh")
          )
        ))
    )
)