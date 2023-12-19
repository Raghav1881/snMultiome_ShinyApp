# Initialize Shiny server
ui <- fluidPage(theme = mtheme,
  titlePanel(
    "snATAC- and snRNA-Seq Atlas of ALS/FTLD Orbitofrontal Cortex"),
    navbarPage(
      NULL,
      theme = mtheme,
      # Create tab panel for UMAP-RNA
      rnaUI("rna"),
      # Create tab panel for UMAP-ATAC
      atacUI("atac"),
      # Create tab panel for Coverage and Violin plots
      coverageUI("coverage"),
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
              width = "100%"
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