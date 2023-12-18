rnaUI <- function(id) {
  tabPanel(
    "UMAP-RNA",
    HTML("Dimensional and Feature Plot of snRNA-seq Data"),
    sidebarLayout(
      sidebarPanel(
        width = 2,
        selectizeInput(NS(id, "genediag1"),
                       "Select gene",
                       choices = NULL),
        selectizeInput(NS(id, "subcatRNA"),
                       "Select identity class",
                       choices = NULL)
      ),
      mainPanel(
        width = 10,
        fluidRow(
          h1("Dimensional plot and Feature Plot"),
          column(
            6,
            plotOutput(NS(id, "dimPlotRNA")),
            downloadButton(NS(id, "dimPlotDownload"),
                           label = "")
          ),
          column(
            6,
            plotOutput(NS(id, "featPlotRNA")),
            downloadButton(NS(id, "featPlotRNADownload"),
                           label = "")
          )
        ),
        fluidRow(
          hr(),
          h1("Feature plot by diagnosis"),
          column(
            12,
            plotOutput(NS(id, "dimPlotRNACtrl")),
            downloadButton(NS(id, "dimPlotRNACtrlDownload"),
                           label = "")
          )
        )
      )
    )
  )
}

rnaServer <- function(id) {
  moduleServer(id, function(input, output, session) {
    updateSelectizeInput(
      session,
      "genediag1",
      choices = geneList,
      server = TRUE,
      selected = "C9orf72"
    )
    updateSelectizeInput(
      session,
      "subcatRNA",
      choices = idents,
      server = TRUE,
      selected = "celltype"
    )
    currentGeneDiag1 <- reactive({
      GetGeneList(req(input$diagchk1),
                  req(input$celltype1))
    })
    output$dimPlotRNA <- renderPlot({
      req(input$genediag1)
      DimPlot(dataset,
              group.by = input$subcatRNA)
    })
    output$dimPlotDownload <- downloadHandler(
      filename = function() {
        paste("DimPlot", input$genediag1, ".png", sep = "")
      },
      content = function(file)  {
        ggsave(file,
               DimPlot(dataset,
                       group.by = input$subcatRNA))
      }
    )
    output$featPlotRNA <- renderPlot({
      req(input$genediag1)
      FeaturePlot(dataset,
                  features = input$genediag1)
    })
    output$featPlotRNADownload <- downloadHandler(
      filename = function() {
        paste("FeatPlot", input$genediag1, ".png", sep = "")
      },
      content = function(file)  {
        ggsave(file,
               FeaturePlot(dataset,
                           features = input$genediag1))
      }
    )
    output$dimPlotRNACtrl <- renderPlot({
      req(input$genediag1)
      FeaturePlot(dataset,
                  features = input$genediag1,
                  split.by = "diagnosis")
    })
    output$dimPlotRNACtrlDownload <- downloadHandler(
      filename = function() {
        paste("FeatPlotDiagnoses", input$genediag1, ".png", sep = "")
      },
      content = function(file) {
        ggsave(
          file,
          FeaturePlot(
            dataset,
            features = input$genediag1,
            split.by = "diagnosis"
          ),
          width = 30
        )
      }
    )
  })
}