library(shiny)
library(Seurat)
library(ggplot2)
library(uwot)
library(plotly)
library(gridExtra)

metadata <- colnames(dataset[[]])
gene_list <- rownames(x = dataset)
diagcell <- vector()
tempdiagcell <- vector()
tempcellinput <- vector()
diagcount <- 0
cellcount <- 0

genConcList <- function(diaginput, cellinput) {
  if (length(diaginput) != length(tempdiagcell) |
  length(cellinput) != length(tempcellinput)) {
    lstAll <- list()
    for (k in length(tempdiagcell)) {
      for (l in length(tempcellinput)) {
        lstAll[k] <- c(paste0(tempdiagcell[k], tempcellinput[l],
                                    collapse = "_"))
      }
    }
  } else {
    validate(
      need(length(diaginput) | length(cellinput) > 0 &
           length(diaginput) == length(tempdiagcell) &
           length(cellinput) == length(tempcellinput),
           "Length of diagnosis and cell types must be equal"))
  }
  return(lstAll)
}

server <- function(input, output, session) {
  # Update dropdown input boxes for all pages
  updateSelectizeInput(session, "select",
                      choices = metadata,
                      server = TRUE)
  updateSelectizeInput(session, "gene_input",
                      choices = gene_list,
                      server = TRUE)
  updateSelectizeInput(session, "genediag",
                      choices = gene_list,
                      server = TRUE)
  updateSelectizeInput(session, "celltype",
                      choices = levels(dataset$celltype),
                      server = TRUE)
gctypes <- reactive(genConcList(input$diagchk, input$celltype))

vlnpt1 <- reactive({
  VlnPlot(dataset,
          features = input$genediag,
          idents = gctypes$temp[1]) &
  theme(axis.title.x = element_blank(),
         legend.position = "None")
})

vlnpt2 <- reactive({
  VlnPlot(dataset,
          features = input$genediag,
          idents = gctypes$temp[2]) &
  theme(axis.title.x = element_blank(),
         legend.position = "None")
})

vlnpt3 <- reactive({
  VlnPlot(dataset,
          features = input$genediag,
          idents = gctypes$temp[3]) &
  theme(axis.title.x = element_blank(),
         legend.position = "None")
})

vlnpt4 <- reactive({
  VlnPlot(dataset,
          features = input$genediag,
          idents = gctypes$temp[4]) &
  theme(axis.title.x = element_blank(),
         legend.position = "None")
})

vlnpt5 <- reactive({
  VlnPlot(dataset,
          features = input$genediag,
          idents = gctypes$temp[5]) &
  theme(axis.title.x = element_blank(),
         legend.position = "None")
})

  # Generate output for features plots
  output$features_graph <- renderPlot({
    validate(
      need(input$feats, "No features selected"))
    VlnPlot(dataset,
      features = input$feats,
      pt.size = 0,
      ncol = 3,
      group.by = "group") &
    theme(axis.title.x = element_blank())
  })

  output$umap_graph <- renderPlot({
    DimPlot(dataset)
  })

  output$categorical_plot <- renderPlot({
    UMAPPlot(dataset, group.by = input$select)
  })

  output$gene_plot <- renderPlot({
    FeaturePlot(dataset, features = input$gene_input)
  })

  # Concatenate strings for use in gene expr plots
  output$vln_gene_plot <- renderPlot({
    ptlist <- list(vlnpt1(), vlnpt2(), vlnpt3(), vlnpt4(), vlnpt5())
    grid.arrange(grobs = ptlist, ncol = length(ptlist))
      })
}