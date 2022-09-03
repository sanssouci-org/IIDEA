#' Open IIDEA
#'
#' @return An object that represents the app. Printing the object or passing it to runApp() will run the app.
#'
#' @importFrom shiny shinyAppDir
#' @export
#'
shinyApp <- function(){
  shinyAppDir("ShinyApps/")
}
