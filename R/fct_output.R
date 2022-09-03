#' Create example data set
#'
#' @details Create example data set from the package sanssouci.data.
#' This dataframe aims at exporting in IIDEA (not to be used by IIDEA as public
#' data set)
#'
#' @return
#' @export
#'
#' @import sanssouci.data
#'
#' @examples
exampleData <- function() {

  #### expression matrix
  data <- list()
  data$matrix <- expr_ALL

  # setProgress(value = 0.4, detail = "Cleaning data ... ")

  categ <- colnames(data$matrix)
  data$categ <- rep(1, length(categ))
  data$categ[which(categ == "NEG")] <- 0

  colnames(data$matrix) <- data$categ


  # setProgress(value = 0.4, detail = "Fc and p.value matrix ... ")

  #### degraded matrix
  dex <- rowWelchTests(data$matrix, data$categ)
  data$degrade <- data.frame("p.value" = dex[["p.value"]], "fc" = dex$estimate)
  rm(dex)

  # gene set matrix
  bioFun <- expr_ALL_GO
  mm <- match(base::rownames(bioFun), base::rownames(data$matrix))
  data$biologicalFunc <- bioFun[mm, ]
  rm(bioFun)
  return(data)
}
