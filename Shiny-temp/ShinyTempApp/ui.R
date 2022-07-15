library(Seurat)
library(shinyWidgets)
library(shiny)
library(bslib)

mtheme <- bs_theme(version = 4, bootswatch = "minty")

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
                                levels(dataset$celltype))),
          mainPanel(
            width = 9, verbatimTextOutput("test"), verbatimTextOutput("test2")))),
      tabPanel("Quality Control")
    )
)