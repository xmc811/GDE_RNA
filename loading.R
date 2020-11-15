
packages <- c("shiny",
              "shinyFiles",
              "shinythemes", 
              "shinyWidgets",
              "tidyverse", "magrittr","DT", 
              "RColorBrewer", "circlize", "ComplexHeatmap", "scales",
              "fgsea",
              "xmcutil")

lapply(packages, require, character.only = TRUE)