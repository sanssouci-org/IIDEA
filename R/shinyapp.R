#' Open IIDEA
#'
#' @return An object that represents the app. Printing the object or passing
#' it to runApp() will run the app.
#'
#' @importFrom shiny shinyAppDir
#' @import shinyjs
#' @import DT
#' @import ggplot2
#' @import htmlwidgets
#' @importFrom plotly renderPlotly
#' @importFrom plotly plotlyOutput
#' @import R.cache
#' @import shinyBS
#' @import shinyjs
#' @importFrom matrixStats rowMaxs rowQuantiles
#' @export
#'
shinyApp <- function() {
  shinyAppDir("ShinyApps/")
}
