library(shiny)
library(Seurat)

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
                      server = TRUE)
  updateSelectizeInput(session, "diagchk",
                      choices = diagnosisList[-1],
                      server = TRUE)

  currentGeneDiag <- reactive({getGeneList(input$diagchk, input$celltype)})
  output$test <- renderPrint({print({currentGeneDiag()[1]})})
  output$test2 <- renderPrint({print({currentGeneDiag()[2]})})
})
