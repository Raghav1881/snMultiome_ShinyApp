# Initialize shiny server
function(input, output, session) {
  rnaServer("rna")
  atacServer("atac")
  updateSelectizeInput(session, "genediag3",
                      choices = geneListATAC,
                      server = TRUE,
                      selected = "C9orf72")
  updateSelectizeInput(session, "diagchk3",
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
    CoveragePlot(ATACdataset,
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
      ggsave(file,
            CoveragePlot(ATACdataset,
                        region = input$genediag3,
                        idents = c(currentGeneDiag3()[[1]],
                                   currentGeneDiag3()[[2]])
                        ) &
                        theme(strip.text.y = element_text(size = 10),
                          title = element_text(size = 15)),
            width = 10,
            height = 15)
    }
  )

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

  output$violin1Download <- downloadHandler(
    filename = function() {
      paste("VlnPlot", input$genediag3, "_", input$celltype3, ".png", sep = "")
    },
    content = function(file) {
      ggsave(file,
            VlnPlot(dataset,
                    features = input$genediag3,
                    idents = currentGeneDiag3()[[1]],
                    pt.size = 0
                  ) &
      theme(axis.title.x = element_blank(),
            legend.position = "none"))
    }
  )

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

  output$violin2Download <- downloadHandler(
    filename = function() {
      paste("VlnPlot", input$genediag3, "_",
            input$celltype3, "2", ".png", sep = "")
    },
    content = function(file) {
      ggsave(file,
            VlnPlot(dataset,
                    features = input$genediag3,
                    idents = currentGeneDiag3()[[2]],
                    pt.size = 0
                  ) &
      theme(axis.title.x = element_blank(),
            legend.position = "none"))
    }
  )

  output$features_graph <- renderPlot({
    validate(
      need(input$feats, "No features selected")
    )
    VlnPlot(dataset,
      features = input$feats,
      pt.size = 0,
      ncol = 2,
      group.by = "celltype") &
    theme(axis.title.x = element_blank())
  })

output$features_graphDownload <- downloadHandler(
    filename = function() {
      paste("FeaturesGraph", ".png", sep = "")
    },
    content = function(file) {
      ggsave(file,
             VlnPlot(dataset,
                    features = input$feats,
                    pt.size = 0,
                    ncol = 2,
                    group.by = "celltype") &
            theme(axis.title.x = element_blank()))
    }
  )

  output$features_graphATAC  <- renderPlot({
    validate(
      need(input$featsATAC, "No features selected")
    )
    FragmentHistogram(ATACdataset,
      group.by = "nucleosome_group") &
    theme(axis.title.x = element_blank())
  })

  output$features_graphATACDownload <- downloadHandler(
    filename = function() {
      paste("FeaturesGraphATAC", ".png", sep = "")
    },
    content = function(file) {
      ggsave(file,
            FragmentHistogram(ATACdataset,
                              group.by = "nucleosome_group") &
            theme(axis.title.x = element_blank()))
    }
  )
}