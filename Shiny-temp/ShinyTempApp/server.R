# Helper function to generate concatenated diagnosis_celltype list
### Vectorized version of GetGeneList haven't tested yet prob won't work
# GetGeneList <- function(diaglist, celllist) {
#   lst_all <- lapply(1:(length(diaglist) + 1), function(k) {
#     paste(if (k == 1) "control" else diaglist[k - 1], celllist, sep = "_")
#   })
#   return(lst_all)
# }

GetGeneList <- function(diaglist, celllist) {
  lst_all <- list()
  temp <- c()
  for (k in 1:(length(diaglist) + 1)) {
    for (l in 1:length(celllist)) {
      if (k == 1) {
        temp[l] <- paste("control", celllist[l], sep = "_")
      } else {
        temp[l] <- paste(diaglist[k - 1], celllist[l], sep = "_")
      }
    }
    lst_all[[k]] <- temp
  }
  return(lst_all)
}

# Initialize shiny server
function(input, output, session) {
  updateSelectizeInput(session, "genediag1",
                      choices = geneList,
                      server = TRUE,
                      selected = "C9orf72")
  updateSelectizeInput(session, "subcatRNA",
                      choices = idents,
                      server = TRUE,
                      selected = "celltype")
  updateSelectizeInput(session, "genediag2",
                      choices = geneListATAC,
                      server = TRUE,
                      selected = "C9orf72")
  updateSelectizeInput(session, "subcatATAC",
                      choices = idents,
                      server = TRUE,
                      selected = "celltype")
  updateSelectizeInput(session, "genediag3",
                      choices = geneListATAC,
                      server = TRUE,
                      selected = "C9orf72")
  updateSelectizeInput(session, "diagchk3",
                      choices = diagnosisList[-1],
                      server = TRUE)

  currentGeneDiag1 <- reactive({
    GetGeneList(input$diagchk1,
                input$celltype1)
  }) %>%
  bindCache(input$diagchk1, input$celltype1)

  currentGeneDiag2 <- reactive({
    GetGeneList(input$diagchk2,
                input$celltype2)
  }) %>%
  bindCache(input$diagchk2, input$celltype2)

  currentGeneDiag3 <- reactive({
    GetGeneList(input$diagchk3,
                input$celltype3)
  }) %>%
  bindCache(input$diagchk3, input$celltype3)

  # Page 1 output begins here
  output$dimPlotRNA <- renderPlot({
    req(input$subcatRNA)
    DimPlot(dataset,
            group.by = input$subcatRNA)
  }) %>%
  bindCache(input$subcatRNA)

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
  }) %>%
  bindCache(input$genediag1)

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
  }) %>%
  bindCache(input$genediag1)

  output$dimPlotRNACtrlDownload <- downloadHandler(
    filename = function() {
      paste("FeatPlotDiagnoses", input$genediag1, ".png", sep = "")
    },
    content = function(file) {
      ggsave(file,
            FeaturePlot(dataset,
                        features = input$genediag1,
                        split.by = "diagnosis"),
            width = 30)
    }
  )

  # Page 2 output begins here
  output$dimPlotATAC <- renderPlot({
    req(input$subcatATAC)
    DimPlot(ATACdataset,
            group.by = input$subcatATAC)
  }) %>%
  bindCache(input$subcatATAC)

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
  }) %>%
  bindCache(input$genediag2)

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
  }) %>%
  bindCache(input$genediag2)

  output$dimPlotATACCtrlDownload <- downloadHandler(
    filename = function() {
      paste("FeatPlotDiagnoses", input$genediag1, ".png", sep = "")
    },
    content = function(file) {
      ggsave(file,
            FeaturePlot(ATACdataset,
                        features = input$genediag2,
                        split.by = "diagnosis"),
              width = 30)
    }
  )

  # Page 3 output begins here
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

  # Page 4 output begins here
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
  }) %>%
  bindCache(input$feats)

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
  }) %>%
  bindCache(input$featsATAC)

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
