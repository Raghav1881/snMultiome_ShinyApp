# Initialize shiny server
function(input, output, session) {
  rnaServer("rna")
  atacServer("atac")
  coverageServer('coverage')
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