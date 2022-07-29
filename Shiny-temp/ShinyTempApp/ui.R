library(Seurat)
library(Signac)
library(shinyWidgets)
library(shiny)
library(bslib)
library(Cairo)

mtheme <- bs_theme(version = 5, bootswatch = "materia")
feature_list <- c("nCount_RNA", "nCount_SCT", "nFeature_RNA",
                  "nFeature_SCT", "percent.mt", "percent.rpl", "percent.rps")
feature_listATAC <- c("nucleosome_group")

ui <- fluidPage(theme = mtheme,
  titlePanel(
    "snATAC- and snRNA-Seq Atlas of ALS/FTLD Orbitofrontal Cortex"),
    navbarPage(
      NULL,
      theme = mtheme,
      tabPanel("UMAP-RNA",
        HTML("Dimensional and Feature Plot of snRNA-seq Data"),
        sidebarLayout(
          sidebarPanel(
            width = 2,
            selectizeInput("genediag1",
                           "Select gene",
                           choices = NULL)),
          mainPanel(
            width = 10,
            fluidRow(
              h1("Dimensional plot and Feature Plot"),
              column(
                6,
                plotOutput("dimPlotRNA")),
              column(
                6,
                plotOutput("featPlotRNA"))),
            fluidRow(
              hr(),
              h1("Feature plot by diagnosis"),
              column(
                12,
                plotOutput("dimPlotRNACtrl")))
        ))),
      tabPanel("UMAP-ATAC",
        HTML("Dimensional and Feature Plot of snATAC-seq Data"),
        sidebarLayout(
          sidebarPanel(
            width = 2,
            selectizeInput("genediag2",
                           "Select gene",
                           choices = NULL)),
          mainPanel(
            width = 10,
            fluidRow(
              h1("Dimensional and Feature Plot of ATAC-seq Data"),
              column(
                6,
                plotOutput("dimPlotATAC")),
              column(
                6,
                plotOutput("featPlotATAC"))),
            fluidRow(
              hr(),
              h1("Feature plot by diagnosis"),
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
            selectizeInput("diagchk3", "Select diagnosis",
                            choices = NULL),
            prettyCheckboxGroup("celltype3", "Select cell types",
                                levels(dataset$celltype),
                                selected = "Oligodendrocytes",
                                icon = icon("check-square"),
                                status = "primary",
                                outline = FALSE,
                                animation = "smooth")),
          mainPanel(
            width = 9,
            h1("Coverage Plot"),
            fluidRow(
              plotOutput("coverage_plot",
                          height = "65vh"),
              hr()),
            fluidRow(
              h1("Violin Plots"),
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
            width = 3,
            prettyCheckboxGroup("feats",
                               "snRNA-seq Features",
                                feature_list,
                                selected = "nCount_RNA",
                                icon = icon("check-square"),
                                status = "primary",
                                outline = FALSE,
                                animation = "smooth"),
            hr(),
            prettyCheckboxGroup("featsATAC",
                                "snATAC-seq Features",
                                feature_listATAC,
                                selected = "nCount_RNA",
                                icon = icon("check-square"),
                                status = "primary",
                                outline = FALSE,
                                animation = "smooth")),
          # Plot features graphs
          mainPanel(
            width = 9,
            fluidRow(
              column(
                6,
                h1("snRNA-seq Feature Plots"),
                plotOutput("features_graph",
                          height = "80vh")),
              column(
                6,
                h1("snATAC-seq Feature Plots"),
                plotOutput("features_graphATAC",
                            height = "80vh"))
              )
          )
        ))
    )
)