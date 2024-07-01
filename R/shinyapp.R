#' Open IIDEA
#'
#' @return An object that represents the app. Printing the object or passing
#' it to runApp() will run the app.
#'
#' @importFrom shiny shinyApp
#' @export
#'
run_IIDEA <- function() {
  # shinyAppDir("inst/ShinyApps/")
  shinyApp(ui = app_ui(), server = app_server)
}
