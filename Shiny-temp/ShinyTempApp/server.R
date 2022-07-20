library(shiny)
library(Seurat)
library(ggplot2)

geneList <- rownames(dataset)
geneListATAC <- rownames(ATACdataset@assays$RNA)
diagnosisList <- levels(dataset@meta.data[["diagnosis"]])

getGeneList <- function(diaglist, celllist) {
  lstAll <- list()
  temp <- c()
  for (k in 1:(length(diaglist) + 1)) {
    for (l in 1:length(celllist)) {
      if (k == 1) {
        temp[l] <- paste("control", celllist[l], sep = "_")
      } else {
        temp[l] <- paste(diaglist[k - 1], celllist[l], sep = "_")
      }
    }
    lstAll[[k]] <- temp
  }
  return(lstAll)
}

shinyServer(function(input, output, session) {
  updateSelectizeInput(session, "genediag1",
                      choices = geneList,
                      server = TRUE,
                      selected = "C9orf72")
  updateSelectizeInput(session, "genediag2",
                      choices = geneListATAC,
                      server = TRUE,
                      selected = "C9orf72")
  updateSelectizeInput(session, "genediag3",
                      choices = geneList,
                      server = TRUE,
                      selected = "C9orf72")
  updateSelectizeInput(session, "diagchk1",
                      choices = diagnosisList[-1],
                      server = TRUE)
  updateSelectizeInput(session, "diagchk2",
                      choices = diagnosisList[-1],
                      server = TRUE)
  updateSelectizeInput(session, "diagchk3",
                      choices = diagnosisList[-1],
                      server = TRUE)

  currentGeneDiag1 <- reactive({getGeneList(input$diagchk1, input$celltype1)})
  currentGeneDiag2 <- reactive({getGeneList(input$diagchk2, input$celltype2)})
  currentGeneDiag3 <- reactive({getGeneList(input$diagchk3, input$celltype3)})

  output$dimPlotRNA <- renderPlot({
    DimPlot(dataset,
            group.by = "celltype")
  })

  output$featPlotRNA <- renderPlot({
    FeaturePlot(dataset,
                features = input$genediag1)
  })

  output$dimPlotRNACtrl <- renderPlot({
    FeaturePlot(dataset,
                features = input$genediag1,
                cells = Cells(input$celltype1),
                split.by = FetchData(dataset,
                                     idents = input$celltype1))
  })

  output$dimPlotRNADiag <- renderPlot({
    FeaturePlot(dataset,
                features = input$genediag1,
                split.by = input$diagchk1)
  })

  output$dimPlotATAC <- renderPlot({
    DimPlot(ATACdataset,
            group.by = "celltype")
  })

  output$featPlotATAC <- renderPlot({
    FeaturePlot(ATACdataset,
                features = input$genediag2)
  })

  output$dimPlotATACCtrl <- renderPlot({
    FeaturePlot(ATACdataset,
                features = input$genediag2,
                cells = Cells(input$celltype2),
                split.by = FetchData(ATACdataset,
                                     idents = input$celltype2))
  })

  output$dimPlotATACDiag <- renderPlot({
    FeaturePlot(ATACdataset,
                features = input$genediag2,
                split.by = input$diagchk2)
  })

  output$coverage_plot <- renderPlot({
    validate(
      need(input$genediag3, "Please select a gene"),
      need(input$diagchk3, "Please select a diagnosis"),
      need(input$celltype3, "Please select the celltypes")
    )
    CoveragePlot(ATACdataset,
                 region = input$genediag3,
                 idents = c(currentGeneDiag3()[[1]],
                            currentGeneDiag3()[[2]]))
  })

  output$violin1 <- renderPlot({
    validate(
      need(input$genediag3, ""),
      need(input$diagchk3, ""),
      need(input$celltype3, "")
    )
    VlnPlot(dataset,
      features = input$genediag3,
      idents = currentGeneDiag3()[[1]],
      pt.size = 0
    ) &
    theme(axis.title.x = element_blank(),
          legend.position = "none")
  })

  output$violin2 <- renderPlot({
    validate(
      need(input$genediag3, ""),
      need(input$diagchk3, ""),
      need(input$celltype3, "")
    )
    VlnPlot(dataset,
      features = input$genediag3,
      idents = currentGeneDiag3()[[2]],
      pt.size = 0
    ) &
    theme(axis.title.x = element_blank(),
          legend.position = "none")
  })

  output$features_graph <- renderPlot({
    validate(
      need(input$feats, "No features selected")
    )
    VlnPlot(dataset,
      features = input$feats,
      pt.size = 0,
      ncol = 3,
      group.by = "group") &
    theme(axis.title.x = element_blank())
  })

})