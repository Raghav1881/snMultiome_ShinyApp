library(shiny)
library(Seurat)
library(ggplot2)
library(uwot)
library(plotly)
library(gridExtra)

metadata <- colnames(dataset[[]])
gene_list <- rownames(x = dataset)
diagcell <- vector()
tempdiagcell <- c()
tempcellinput <- c()
diagcount <- 0
cellcount <- 0

genConcList <- function(diaginput, cellinput) {
  if (length(diaginput) != length(tempdiagcell) |
  length(cellinput) != length(tempcellinput)) {
    lstAll <- list()
    tempdiagcell <- input$diagchk
    tempcellinput <- input$celltype
    for (k in 1:length(tempdiagcell)) {
      for (l in 1:length(tempcellinput)) {
        temp[l] <- c(paste(tempdiagcell[k], tempcellinput[l],
                             sep = "_"))
      }
      lstAll[[k]] <- temp
    }
    return(lstAll)
  } else {
    validate(
      need(length(diaginput) == length(tempdiagcell) &
           length(cellinput) == length(tempcellinput),
           "Length of diagnosis and cell types must be equal"))
  }
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
  output$vln_plot1 <- renderPlot({
    gctypes <<- genConcList(input$diagchk, input$celltype)
    VlnPlot(dataset,
            features = input$genediag,
            idents = gctypes[[1]]) &
    theme(axis.title.x = element_blank(),
          legend.position = "None")
})
  output$vln_plot2 <- renderPlot({
    VlnPlot(dataset,
            features = input$genediag,
            idents = gctypes[[2]]) &
    theme(axis.title.x = element_blank(),
          legend.position = "None")
})
  output$vln_plot3 <- renderPlot({
    VlnPlot(dataset,
            features = input$genediag,
            idents = gctypes[[3]]) &
    theme(axis.title.x = element_blank(),
          legend.position = "None")
})
  output$vln_plot4 <- renderPlot({
    VlnPlot(dataset,
            features = input$genediag,
            idents = gctypes[[4]]) &
    theme(axis.title.x = element_blank(),
          legend.position = "None")
})
  output$vln_plot5 <- renderPlot({
    VlnPlot(dataset,
            features = input$genediag,
            idents = gctypes[[5]]) &
    theme(axis.title.x = element_blank(),
          legend.position = "None")
})
}