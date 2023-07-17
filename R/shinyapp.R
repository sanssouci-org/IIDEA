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
#' @import plotly
#' @import R.cache
#' @import shinyBS
#' @import shinyjs
#' @export
#'
shinyApp <- function() {
  shinyAppDir("ShinyApps/")
}
