library("shiny")
library("shinyjs")
library("sanssouci")
library("sanssouci.data")
library("plotly")
stopifnot(packageVersion("sanssouci.data") >= "0.5.0")
library("ggplot2")
library("dplyr")
library("htmlwidgets")
library("DT")
library("shinyBS")
library("stringr")
library("R.cache")
library("IIDEA")



# data(expr_ALL, package = "sanssouci.data", envir = environment())
# data(expr_ALL_GO, package = "sanssouci.data", envir = environment())
# data(RNAseq_blca, package = "sanssouci.data", envir = environment())
# data(RNAseq_blca_GO, package = "sanssouci.data", envir = environment())


IIDEA::app_ui()
