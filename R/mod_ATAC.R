atacUI <- function(id) {
  tabPanel(
    "UMAP-ATAC",
    HTML("Dimensional and Feature Plot of snATAC-seq Data"),
    sidebarLayout(
      sidebarPanel(
        width = 2,
        selectizeInput("genediag2",
                       "Select gene",
                       choices = NULL),
        selectizeInput("subcatATAC",
                       "Select identity class",
                       choices = NULL)
      ),
      mainPanel(
        width = 10,
        fluidRow(
          h1("Dimensional and Feature Plot of ATAC-seq Data"),
          column(
            6,
            plotOutput("dimPlotATAC"),
            downloadButton("dimPlotATACDownload",
                           label = "")
          ),
          column(
            6,
            plotOutput("featPlotATAC"),
            downloadButton("featPlotATACDownload",
                           label = "")
          )
        ),
        fluidRow(
          hr(),
          h1("Feature plot by diagnosis"),
          column(
            12,
            plotOutput("dimPlotATACCtrl"),
            downloadButton("dimPlotATACCtrlDownload",
                           label = "")
          )
        )
      )
    )
  )
}

atacServer <- function(id) {
  moduleServer(id, function(input, output, session) {
    updateSelectizeInput(
      session,
      "genediag2",
      choices = geneListATAC,
      server = TRUE,
      selected = "C9orf72"
    )
    updateSelectizeInput(
      session,
      "subcatATAC",
      choices = idents,
      server = TRUE,
      selected = "celltype"
    )
    currentGeneDiag2 <- reactive({
      GetGeneList(input$diagchk2,
                  input$celltype2)
    })
    output$dimPlotATAC <- renderPlot({
      req(input$genediag2)
      DimPlot(ATACdataset,
              group.by = input$subcatATAC)
    })
    
    output$dimPlotATACDownload <- downloadHandler(
      filename = function() {
        paste("DimPlotATAC", input$genediag1, ".png", sep = "")
      },
      content = function(file) {
        ggsave(file,
               DimPlot(ATACdataset,
                       group.by = input$subcatATAC))
      }
    )
    output$featPlotATAC <- renderPlot({
      FeaturePlot(ATACdataset,
                  features = input$genediag2)
    })
    output$featPlotATACDownload <- downloadHandler(
      filename = function() {
        paste("FeatPlotATAC", input$genediag1, ".png", sep = "")
      },
      content = function(file) {
        ggsave(file,
               FeaturePlot(ATACdataset,
                           features = input$genediag2))
      }
    )
    output$dimPlotATACCtrl <- renderPlot({
      FeaturePlot(ATACdataset,
                  features = input$genediag2,
                  split.by = "diagnosis")
    })
    output$dimPlotATACCtrlDownload <- downloadHandler(
      filename = function() {
        paste("FeatPlotDiagnoses", input$genediag1, ".png", sep = "")
      },
      content = function(file) {
        ggsave(
          file,
          FeaturePlot(
            ATACdataset,
            features = input$genediag2,
            split.by = "diagnosis"
          ),
          width = 30
        )
      }
    )
  })
}