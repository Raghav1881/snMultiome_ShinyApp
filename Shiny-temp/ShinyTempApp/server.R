library(shiny)
library(Seurat)
library(ggplot2)

geneList <- rownames(dataset)
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
  updateSelectizeInput(session, "diagchk1",
                      choices = diagnosisList[-1],
                      server = TRUE)
  updateSelectizeInput(session, "genediag2",
                      choices = geneList,
                      server = TRUE,
                      selected = "C9orf72")
  updateSelectizeInput(session, "diagchk2",
                      choices = diagnosisList[-1],
                      server = TRUE)

  currentGeneDiag1 <- reactive({getGeneList(input$diagchk1, input$celltype1)})
  currentGeneDiag2 <- reactive({getGeneList(input$diagchk2, input$celltype2)})
  output$test <- renderPrint({print({currentGeneDiag()[1]})})
  output$test2 <- renderPrint({print({currentGeneDiag()[2]})})

  output$dimPlotRNA <- renderPlot({
    DimPlot(dataset, group.by = "celltype")
  })

  output$featPlotRNA <- renderPlot({
    FeaturePlot(dataset, features = input$genediag1)
  })

  output$dimPlotRNACtrl <- renderPlot({
    DimPlot(dataset, group.by = currentGeneDiag1()[1])
  })

  output$dimPlotRNADiag <- renderPlot({
    DimPlot(dataset, split.by = currentGeneDiag1()[2])
  })

  output$violin1 <- renderPlot({
    validate(
      need(input$genediag2, "Please select a gene"),
      need(input$diagchk2, "Please select a diagnosis"),
      need(input$celltype2, "Please select the celltypes")
    )
    VlnPlot(dataset,
      features = input$genediag2,
      idents = currentGeneDiag2()[[1]],
      pt.size = 0
    ) &
    theme(axis.title.x = element_blank(),
          legend.position = "none")
  })

  output$violin2 <- renderPlot({
    validate(
      need(input$genediag2, ""),
      need(input$diagchk2, ""),
      need(input$celltype2, "")
    )
    VlnPlot(dataset,
      features = input$genediag2,
      idents = currentGeneDiag2()[[2]],
      pt.size = 0
    ) &
    theme(axis.title.x = element_blank(),
          legend.position = "none")
  })

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

})
