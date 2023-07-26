#' Post hoc confidence bounds on the true/false positives from An object
#' of class 'SansSouci'
#'
#' Imported from sanssouci package.
#'
#' @param object An object of class 'SansSouci'
#' @param S  A subset of indices
#' @param what A character vector, the names of the post hoc bounds to be
#'   computed, among:
#'
#' - FP: Upper bound on the number of false positives in the 'x' most
#' significant items
#' - TP: Lower bound on the number of true positives in the 'x' most
#' significant items
#' - FDP: Upper bound on the proportion of false positives in the 'x' most
#'  significant items
#' - TP: Lower bound on the proportion of true positives in the 'x' most
#' significant items
#'
#' Defaults to `c("TP", "FDP")`
#' @param all A logical value: should the bounds for all ordered subsets of `S`
#'  be returned? If `FALSE` (the default), only the bound for `S` is returned
#' @param ... Not used
#'
#' @return If `all` is `FALSE` (the default), only the value of the bound is
#'  returned. Otherwise, a `data.frame` is return, with |S| rows and 4 columns:
#' - x: Number of most significant items selected
#' - label: Label for the procedure, typically of the form 'family(param)'
#' - bound: Value of the post hoc bound
#' - stat: Type of post hoc bound, as specified by argument `bound`.
#'
#' @importFrom stats predict
predict2 <- function(object, S = seq_len(nHyp(object)),
                     what = c("TP", "FDP"), all = FALSE, ...) {
  p.values <- pValues(object)
  thr <- thresholds(object)
  lab <- label(object)
  if (max(S) > nHyp(object)) {
    stop("'S' is not a subset of hypotheses")
  }
  bounds <- posthoc_bound2(p.values, S = S, thr = thr, lab = lab,
                           what = what, all = all)
  if (!all) {
    bounds <- bounds[, "bound"]
    if (length(bounds) > 1) {
      names(bounds) <- what
    }
  }
  return(bounds)
}

#' Post hoc confidence bounds on the true/false positives from pvalues
#'
#' Imported from sanssouci package.
#'
#' @param p.values a matrix containing m x B pvalues
#' @param S  A subset of indices
#' @param thr a vector of calibration thresholds
#' @param lab charactor
#' @param what A character vector, the names of the post hoc bounds to be
#'   computed, among:
#'
#' - FP: Upper bound on the number of false positives in the 'x' most
#' significant items
#' - TP: Lower bound on the number of true positives in the 'x' most
#'  significant items
#' - FDP: Upper bound on the proportion of false positives in the 'x' most
#' significant items
#' - TP: Lower bound on the proportion of true positives in the 'x' most
#'  significant items
#'
#' Defaults to `c("TP", "FDP")`
#' @param all A logical value: should the bounds for all ordered subsets of `S`
#'  be returned? If `FALSE` (the default), only the bound for `S` is returned
#'
#' @return If `all` is `FALSE` (the default), only the value of the bound is
#'  returned. Otherwise, a `data.frame` is return, with |S| rows and 4 columns:
#' - x: Number of most significant items selected
#' - label: Label for the procedure, typically of the form 'family(param)'
#' - bound: Value of the post hoc bound
#' - stat: Type of post hoc bound, as specified by argument `bound`.
#'
#' @export
#' @importFrom stats predict
#' @import sanssouci
posthoc_bound2 <- function(p.values, S = seq_along(p.values),
                           thr = NULL, lab = NULL,
                           what = c("TP", "FDP"), all = FALSE) {
  if (is.null(thr)) {
    stop("Argument 'thr' must be non NULL")
  }
  if (is.null(lab)) {
    lab <- seq(1, length(p.values))
  }
  s <- length(S)
  idxs <- seq_len(s)
  max_FP <- rep(NA_integer_, s)
  pS <- p.values[S]
  o <- order(pS)
  sorted_p <- pS[o]
  if (length(thr) == length(p.values) && all(thr %in% c(0, 1))) {
    max_FP <- cumsum(thr[o] == 0)
  } else {
    if (s <= 5 * sqrt(length(thr))) {
      max_FP <- maxFP(sorted_p, thr)
      idxs <- length(idxs)
    } else {
      max_FP <- curveMaxFP(sorted_p, thr)
    }
  }
  bounds <- formatBounds(max_FP,
    idxs = idxs, lab = lab, what = what,
    all = all
  )
  bounds
}


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
boundGroup2 <- function(object) {
  table <- data.frame("Name" = c(), "# genes" = c(), "TP&ge;" = c(),
                      "FDP&le;" = c(), check.names = FALSE)
  bioFun <- object$input$biologicalFunc
  print(paste("boundGroup2 class is ", class(bioFun)))
  if (class(bioFun)[1] == "list") {
    nameFunctions <- names(bioFun)
    for (func in nameFunctions) {
      incProgress(1 / length(nameFunctions))
      name <- bioFun[[func]]
      ids <- which(is.element(rownames(object$input$Y), name))
      if (length(ids) > 5) {
        bounds <- predict2(object, S = ids)
        table <- rbind(table, data.frame(
          "Name" = addUrlLink(func),
          "# genes" = length(ids),
          "TP&ge;" = as.integer(bounds["TP"]),
          "FDP&le;" = bounds["FDP"],
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
        bounds <- predict2(object, S = ids)
        table <- rbind(table, data.frame(
          "Name" = addUrlLink(func),
          "# genes" = length(ids),
          "TP&ge;" = as.integer(bounds["TP"]),
          "FDP&le;" = bounds["FDP"],
          check.names = FALSE, row.names = NULL
        ))
      }
    }
  }
  table <- table[order(table["FDP&le;"]), ]
  print("end boundGroup2")
  print("###################")
  return(table)
}
