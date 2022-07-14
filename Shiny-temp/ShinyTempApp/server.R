library(shiny)
library(Seurat)
geneList <- rownames(dataset)

# Define server logic required to draw a histogram
shinyServer(function(input, output, session) {
  updateSelectizeInput(session, "genediag",
                       choices = geneList,
                       server = TRUE)

})
