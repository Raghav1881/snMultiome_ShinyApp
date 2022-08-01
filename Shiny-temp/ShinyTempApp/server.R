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
                      choices = geneListATAC,
                      server = TRUE,
                      selected = "C9orf72")
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

  output$dimPlotDownload <- downloadHandler(
    filename = function() {
      paste("DimPlot", input$genediag1, ".png", sep = "")
    },
    content = function(file)  {
      ggsave(file,
            DimPlot(dataset,
                    group.by = "celltype"))
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
    content = function(file)  {
      ggsave(file,
            FeaturePlot(dataset,
                        features = input$genediag1,
                        split.by = "diagnosis"))
    }
  )

  output$dimPlotATAC <- renderPlot({
    DimPlot(ATACdataset,
            group.by = "celltype")
  })

  output$dimPlotATACDownload <- downloadHandler(
    filename = function() {
      paste("DimPlotATAC", input$genediag1, ".png", sep = "")
    },
    content = function(file)  {
      ggsave(file,
            DimPlot(ATACdataset,
                    group.by = "celltype"))
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
    content = function(file)  {
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
    content = function(file)  {
      ggsave(file,
            FeaturePlot(ATACdataset,
                        features = input$genediag2,
                        split.by = "diagnosis"))
    }
  )

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
      paste("FeatPlotDiagnoses", input$genediag1, ".png", sep = "")
    },
    content = function(file)  {
      ggsave(file,
            CoveragePlot(ATACdataset,
                        region = input$genediag3,
                        idents = c(currentGeneDiag3()[[1]],
                                   currentGeneDiag3()[[2]])
                        ) &
                        theme(strip.text.y = element_text(size = 10),
                          title = element_text(size = 15)))
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
    content = function(file)  {
      ggsave(file,
            VlnPlot(dataset,
                    features = input$genediag3,
                    idents = currentGeneDiag3()[[1]],
                    pt.size = 0
                  ) &
      theme(axis.title.x = element_blank(),
            legend.position = "none"))
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

output$violin2Download <- downloadHandler(
    filename = function() {
      paste("VlnPlot", input$genediag3, "_", input$celltype3, "2", ".png", sep = "")
    },
    content = function(file)  {
      ggsave(file,
            VlnPlot(dataset,
                    features = input$genediag3,
                    idents = currentGeneDiag3()[[2]],
                    pt.size = 0
                  ) &
      theme(axis.title.x = element_blank(),
            legend.position = "none"))
  })

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
    content = function(file)  {
      ggsave(file,
             VlnPlot(dataset,
                    features = input$feats,
                    pt.size = 0,
                    ncol = 2,
                    group.by = "celltype") &
            theme(axis.title.x = element_blank()))
  })

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
    content = function(file)  {
      ggsave(file,
            FragmentHistogram(ATACdataset,
                              group.by = "nucleosome_group") &
            theme(axis.title.x = element_blank()))
  })
})