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
  updateSelectizeInput(session, "genediag",
                      choices = geneList,
                      server = TRUE,
                      selected = "C9orf72")
  updateSelectizeInput(session, "diagchk",
                      choices = diagnosisList[-1],
                      server = TRUE)

  currentGeneDiag <- reactive({getGeneList(input$diagchk, input$celltype)})
  output$test <- renderPrint({print({currentGeneDiag()[1]})})
  output$test2 <- renderPrint({print({currentGeneDiag()[2]})})

  output$violin1 <- renderPlot({
    validate(
      need(input$genediag, "Please select a gene"),
      need(input$diagchk, "Please select a diagnosis"),
      need(input$celltype, "Please select the celltypes")
    )
    VlnPlot(dataset,
      features = input$genediag,
      idents = currentGeneDiag()[[1]],
      pt.size = 0
    ) &
    theme(axis.title.x = element_blank())
  })

  output$violin2 <- renderPlot({
    validate(
      need(input$genediag, ""),
      need(input$diagchk, ""),
      need(input$celltype, "")
    )
    VlnPlot(dataset,
      features = input$genediag,
      idents = currentGeneDiag()[[2]],
      pt.size = 0
    ) &
    theme(axis.title.x = element_blank())
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
