library(Seurat)
library(shiny)
library(stringr)

diagnosisList <- levels(dataset@meta.data[["diagnosis"]])

ui <- fluidPage(
  titlePanel(
    "scRNA-seq Brain ALS Data"),
    navbarPage(
      NULL,
      tabPanel("snRNA-seq",
         HTML("Diagnosis and Gene Expression Analysis"),
         sidebarLayout(
           sidebarPanel(
             width = 3,
             selectizeInput("genediag", "Select gene", choices = NULL),
             pickerInput("diagchk",
                         label = "Select diagnosis",
                        choices = diagnosisList[-1]),
             checkboxGroupInput("celltype", "Select cell types",
                                levels(dataset$celltype))),
           # Generate violin plot output based on gene, diagnosis, and cell types
           mainPanel(
             width = 9, plotOutput("vln_gene_plot"))))
    )
)

