#' IIDEA post-hoc bound table on gene set.
#'
#' @param object An object of class 'SansSouci'
#'
#' @return a data.frame containing:
#'
#' - the name of the gene set
#' - the number of genes in the gene set
#' - the estimate TP bound
#' - the estimate FDP bound
#'
#' @export
#' @importFrom shiny incProgress
#' @import sanssouci
boundGroup2 <- function(object) {
  table <- data.frame(
    "Name" = c(), "# genes" = c(), "TP\u2265" = c(),
    "FDP\u2264" = c(), check.names = FALSE
  )
  bioFun <- object$input$biologicalFunc
  if (class(bioFun)[1] == "list") {
    nameFunctions <- names(bioFun)
    for (func in nameFunctions) {
      incProgress(1 / length(nameFunctions))
      name <- bioFun[[func]]
      ids <- which(is.element(rownames(object$input$Y), name))
      if (length(ids) > 5) {
        bounds <- predict(object, S = ids)
        table <- rbind(table, data.frame(
          "Name" = addUrlLink(func),
          "# genes" = length(ids),
          "TP\u2265" = as.integer(bounds["TP"]),
          "FDP\u2264" = bounds["FDP"],
          check.names = FALSE, row.names = NULL
        ))
      }
    }
  } else {
    nameFunctions <- colnames(bioFun)
    for (func in nameFunctions) {
      # incProgress(1 / length(nameFunctions))
      ids <- which(bioFun[, func] == 1)
      if (length(ids) > 1) {
        bounds <- predict(object, S = ids)
        table <- rbind(table, data.frame(
          "Name" = addUrlLink(func),
          "# genes" = length(ids),
          "TP\u2265" = as.integer(bounds["TP"]),
          "FDP\u2264" = bounds["FDP"],
          check.names = FALSE, row.names = NULL
        ))
      }
    }
  }
  table <- table[order(table["FDP\u2264"]), ]
  return(table)
}
