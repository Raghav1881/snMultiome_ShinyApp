# global.R file used to load common datasets shared across all sessions
setwd("~/Documents/GitHub/LMP-RShiny-App/Shiny-temp/ShinyTempApp/data/")
if (!require("pacman")) {
  install.packages("pacman")
}

library(pacman)

p_load("shiny", "Seurat", "Signac", "shinyWidgets",
       "bslib", "ggplot2")

dataset <- readRDS("30k_snRNA.RDS")
ATACdataset <- readRDS("9k_snATAC.RDS")
geneList <- rownames(dataset)
geneListATAC <- rownames(ATACdataset@assays$RNA)
diagnosisList <- levels(dataset@meta.data[["diagnosis"]])
mtheme <- bs_theme(version = 5, bootswatch = "materia")
feature_list <- c("nCount_RNA", "nCount_SCT", "nFeature_RNA",
                  "nFeature_SCT", "percent.mt", "percent.rpl", "percent.rps")
feature_listATAC <- c("nucleosome_group")
idents <- c("sample", "sex", "celltype", "diagnosis",
            "cellsubtype")