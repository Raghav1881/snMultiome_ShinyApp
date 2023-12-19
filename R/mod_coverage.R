coverageUI <- function(id) {
  tabPanel(
    "Coverage/Violin Plots",
    HTML("Diagnosis and Gene Expression Analysis"),
    sidebarLayout(
      sidebarPanel(
        width = 3,
        selectizeInput(NS(id, "genediag3"), "Select gene",
                       choices = NULL),
        selectizeInput(NS(id, "diagchk3"), "Select diagnosis",
                       choices = NULL),
        prettyCheckboxGroup(
          NS(id, "celltype3"),
          "Select cell types",
          levels(dataset$celltype),
          selected = "Oligodendrocytes",
          icon = icon("check-square"),
          status = "primary",
          outline = FALSE,
          animation = "smooth"
        )
      ),
      mainPanel(
        # Plot coverage plot for all diagnoses/celltypes
        width = 9,
        h1("Coverage Plot"),
        fluidRow(plotOutput(NS(id, "coverage_plot"),
                            height = "65vh")),
        fluidRow(column(
          6,
          downloadButton("coverage_plotDownload",
                         label = "")
        ),
        hr()),
        # Plot violin plots for each diagnosis
        fluidRow(
          h1("Violin Plots"),
          column(
            6,
            # Violin plot for control_celltype
            plotOutput(NS(id, "violin1"),
                       height = "40vh"),
            downloadButton("violin1Download",
                           label = "")
          ),
          column(
            6,
            # Violin plot for diagnosis_celltype
            plotOutput(NS(id, "violin2"),
                       height = "40vh"),
            downloadButton("violin2Download",
                           label = "")
          )
        )
      )
    )
  )
}
coverageServer <- function(id) {
  moduleServer(id, function(input, output, session) {
    updateSelectizeInput(
      session,
      "genediag3",
      choices = geneListATAC,
      server = TRUE,
      selected = "C9orf72"
    )
    updateSelectizeInput(session,
                         "diagchk3",
                         choices = diagnosisList[-1],
                         server = TRUE)
    currentGeneDiag3 <- reactive({
      GetGeneList(input$diagchk3,
                  input$celltype3)
    })
    
    output$coverage_plot <- renderPlot({
      validate(
        need(input$genediag3, "Please select a gene"),
        need(input$diagchk3, "Please select a diagnosis"),
        need(input$celltype3, "Please select the celltypes")
      )
      CoveragePlot(
        ATACdataset,
        region = input$genediag3,
        idents = c(currentGeneDiag3()[[1]],
                   currentGeneDiag3()[[2]])
      ) &
        theme(strip.text.y = element_text(size = 10),
              title = element_text(size = 15))
    })
    
    output$coverage_plotDownload <- downloadHandler(
      filename = function() {
        paste("CoveragePlot", input$genediag1, ".png", sep = "")
      },
      content = function(file) {
        ggsave(
          file,
          CoveragePlot(
            ATACdataset,
            region = input$genediag3,
            idents = c(currentGeneDiag3()[[1]],
                       currentGeneDiag3()[[2]])
          ) &
            theme(
              strip.text.y = element_text(size = 10),
              title = element_text(size = 15)
            ),
          width = 10,
          height = 15
        )
      }
    )
    
    output$violin1 <- renderPlot({
      validate(
        need(input$genediag3, ""),
        need(input$diagchk3, ""),
        need(input$celltype3, "")
      )
      VlnPlot(
        dataset,
        features = input$genediag3,
        idents = currentGeneDiag3()[[1]],
        pt.size = 0
      ) &
        theme(axis.title.x = element_blank(),
              legend.position = "none")
    })
    
    output$violin1Download <- downloadHandler(
      filename = function() {
        paste("VlnPlot",
              input$genediag3,
              "_",
              input$celltype3,
              ".png",
              sep = "")
      },
      content = function(file) {
        ggsave(
          file,
          VlnPlot(
            dataset,
            features = input$genediag3,
            idents = currentGeneDiag3()[[1]],
            pt.size = 0
          ) &
            theme(axis.title.x = element_blank(),
                  legend.position = "none")
        )
      }
    )
    
    output$violin2 <- renderPlot({
      validate(
        need(input$genediag3, ""),
        need(input$diagchk3, ""),
        need(input$celltype3, "")
      )
      VlnPlot(
        dataset,
        features = input$genediag3,
        idents = currentGeneDiag3()[[2]],
        pt.size = 0
      ) &
        theme(axis.title.x = element_blank(),
              legend.position = "none")
    })
    
    output$violin2Download <- downloadHandler(
      filename = function() {
        paste("VlnPlot",
              input$genediag3,
              "_",
              input$celltype3,
              "2",
              ".png",
              sep = "")
      },
      content = function(file) {
        ggsave(
          file,
          VlnPlot(
            dataset,
            features = input$genediag3,
            idents = currentGeneDiag3()[[2]],
            pt.size = 0
          ) &
            theme(axis.title.x = element_blank(),
                  legend.position = "none")
        )
      }
    )
  })
}