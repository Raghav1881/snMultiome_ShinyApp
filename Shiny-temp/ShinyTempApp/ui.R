# Initialize Shiny server
ui <- fluidPage(theme = mtheme,
  titlePanel(
    "snATAC- and snRNA-Seq Atlas of ALS/FTLD Orbitofrontal Cortex"),
    navbarPage(
      NULL,
      theme = mtheme,
      # Create tab panel for UMAP-RNA
      tabPanel("UMAP-RNA",
        HTML("Dimensional and Feature Plot of snRNA-seq Data"),
        sidebarLayout(
          sidebarPanel(
            width = 2,
            selectizeInput("genediag1",
                           "Select gene",
                           choices = NULL),
            selectizeInput("subcatRNA",
                           "Select identity class",
                           choices = NULL)),
          mainPanel(
            width = 10,
            fluidRow(
              h1("Dimensional plot and Feature Plot"),
              column(
                6,
                plotOutput("dimPlotRNA"),
                downloadButton("dimPlotDownload",
                                label = "")),
              column(
                6,
                plotOutput("featPlotRNA"),
                downloadButton("featPlotRNADownload",
                                label = ""))),
            fluidRow(
              hr(),
              h1("Feature plot by diagnosis"),
              column(
                12,
                plotOutput("dimPlotRNACtrl"),
                downloadButton("dimPlotRNACtrlDownload",
                                label = "")))
        ))),
      # Create tab panel for UMAP-ATAC
      tabPanel("UMAP-ATAC",
        HTML("Dimensional and Feature Plot of snATAC-seq Data"),
        sidebarLayout(
          sidebarPanel(
            width = 2,
            selectizeInput("genediag2",
                           "Select gene",
                           choices = NULL),
            selectizeInput("subcatATAC",
                           "Select identity class",
                           choices = NULL)),
          mainPanel(
            width = 10,
            fluidRow(
              h1("Dimensional and Feature Plot of ATAC-seq Data"),
              column(
                6,
                plotOutput("dimPlotATAC"),
                downloadButton("dimPlotATACDownload",
                                label = "")),
              column(
                6,
                plotOutput("featPlotATAC"),
                downloadButton("featPlotATACDownload",
                                label = ""))),
            fluidRow(
              hr(),
              h1("Feature plot by diagnosis"),
              column(
                12,
                plotOutput("dimPlotATACCtrl"),
                downloadButton("dimPlotATACCtrlDownload",
                                label = "")))

        ))),
      # Create tab panel for Coverage and Violin plots
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
            # Plot coverage plot for all diagnoses/celltypes
            width = 9,
            h1("Coverage Plot"),
            fluidRow(
              plotOutput("coverage_plot",
                          height = "65vh")),
            fluidRow(
              column(
                6,
                downloadButton("coverage_plotDownload",
                              label = "")),
              hr()),
            # Plot violin plots for each diagnosis
            fluidRow(
              h1("Violin Plots - RNA"),
              column(
                6,
                # Violin plot for control_celltype
                plotOutput("violin1",
                            height = "40vh"),
                downloadButton("violin1Download",
                                label = "")),
              column(
                6,
                # Violin plot for diagnosis_celltype
                plotOutput("violin2",
                            height = "40vh"),
                downloadButton("violin2Download",
                                label = "")))
            )
        )),
      # Create tab panel for Quality Control
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
                                selected = "nucleosome_group",
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
                          height = "80vh"),
                downloadButton("features_graphDownload",
                                label = "")),
              column(
                6,
                h1("snATAC-seq Feature Plots"),
                plotOutput("features_graphATAC",
                            height = "80vh"),
                downloadButton("features_graphATACDownload",
                                label = "")))
          )
        )),
      tabPanel("References",
        fluidRow(
          column(6,
            h1("Acknowledgements"),
            tags$img(
              src = "Affiliations.png",
              width = "80%"
            )
          ),
          column(
            6,
            h1("Packages"),
            tagList(a("Seurat",
                      href = "https://satijalab.org/seurat/index.html")),
            p("Hao Y, Hao S, Andersen-Nissen E, Mauck III WM, Zheng S,
              Butler A, Lee MJ, Wilk AJ, Darby C, Zager M, Hoffman P.
              Integrated analysis of multimodal single-cell data.
              Cell. 2021 Jun 24;184(13):3573-87."),
            tagList(a("Signac",
                      href = "https://satijalab.org/signac/index.html")),
            p("Stuart et al.
              Single-cell chromatin state analysis with Signac.
              Nature Methods. 2021"),
            tagList(a("Shiny",
                      href = "https://CRAN.R-project.org/package=shiny")),
            p("Chang W, Cheng J, Allaire JJ, Sievery C, Schloerke B, Xie Y,
              Allen J, Mcpherson J, Dipert A, Borges B.
              shiny: Web Application Framework for R. 2022."),
            tagList(a("ggplot2",
                      href = "https:/ggplot2.tidyverse.org")),
            p("Wickham H. ggplot2: Elegant Graphics for Data Analysis.
              Springer-Verlag New York. 2016."),
            tagList(a("bslib",
                      href = "https://CRAN.R-project.org/package=bslib")),
            p("Sievery C, Cheng J.
              bslib: Custom Bootstrap Sass Themes for shiny and rmarkdown.
              2022.")
          )
        )
      )
    )
)