library(shiny)
library(Seurat)
library(ggplot2)
library(uwot)
library(plotly)

metadata <- colnames(dataset[[]])
gene_list <- rownames(x = dataset)

genConcList <- function(diaginput, cellinput) {
  if (length(diaginput) == length(cellinput)) {
    for (k in diaginput) {
      tmp <- c(diaginput, cellinput)
      diagcell[k] <- paste0(tmp, collapse = "_")
    }
  } else {
      validate(
        need(length(diaginput) == length(cellinput),
            "Length of diagnosis and cell types must be equal"))
  }
  return(diagcell)
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
    gctypes <- genConcList(input$diagchk, input$celltype)
    VlnPlot(dataset,
            features = input$genediag,
            idents = gctypes) &
    theme(axis.title.x = element_blank(),
          legend.position = "None")
  })
}