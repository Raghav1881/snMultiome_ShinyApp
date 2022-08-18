setwd("~/Documents/GitHub/LMP-RShiny-App/Shiny-temp/ShinyTempApp/data/")
dataset <- readRDS("30k_snRNA.RDS")
ATACdataset <- readRDS("9k_snATAC.RDS")
p_load("shiny", "Seurat", "Signac", "shinyWidgets",
       "bslib", "ggplot2")

if (!require("pacman")) {
  install.packages("pacman")
}
