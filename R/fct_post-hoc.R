#' Upper bound for the number of false discoveries among most significant items
#'
#' Imported from sanssouci package
#'
#' @param p.values A vector containing s p-values
#' @param thr A vector of \eqn{K} JER-controlling thresholds

#' @return A vector of size \eqn{s} giving an joint upper confidence bound on
#'   the number of false discoveries among the \eqn{k} most significant items
#'   for all \eqn{i \in \{1\ldots s\}}.

#' @author Gilles Blanchard, Nicolas Enjalbert-Courrech, Pierre Neuvial and
#'   Etienne Roquain
#' @details The time and space complexity of this function is O(s), which is
#'   optimal since s is the length of the returned vector.
curveMaxFP <- function(p.values, thr) {
  s <- length(p.values)
  if (s == 0) {
    return(numeric(0L))
  }
  p.values <- sort(p.values)
  thr <- sort(thr)

  kMax <- length(thr)
  if (s < kMax) { # truncate thr to first 's' values
    seqK <- seq(from = 1, to = s, by = 1)
    thr <- thr[seqK]
  } else { # complete 'thr' to length 's' with its last value
    thr <- c(thr, rep(thr[kMax], s - kMax))
  }
  ## sanity checks
  stopifnot(length(thr) == s)
  rm(kMax)

  K <- rep(s, s) ## K[i] = number of k/ T[i] <= s[k]
  Z <- rep(s, s) ## Z[k] = number of i/ T[i] >  s[k] = cardinal of R_k
  ## 'K' and 'Z' are initialized to their largest possible value (both 's')
  kk <- 1
  ii <- 1
  while ((kk <= s) && (ii <= s)) {
    if (thr[kk] > p.values[ii]) {
      K[ii] <- kk - 1
      ii <- ii + 1
    } else {
      Z[kk] <- ii - 1
      kk <- kk + 1
    }
  }
  Vbar <- numeric(s)
  ww <- which(K > 0)
  A <- Z - (1:s) + 1
  cA <- cummax(A)[K[ww]] # cA[i] = max_{k<K[i]} A[k]
  Vbar[ww] <- pmin(ww - cA, K[ww])
  Vbar
}

#' Post hoc confidence bounds on the true/false positives from An object of class 'SansSouci'
#'
#' Imported from sanssouci package.
#'
#' @param object An object of class 'SansSouci'
#' @param S  A subset of indices
#' @param what A character vector, the names of the post hoc bounds to be
#'   computed, among:
#'
#' - FP: Upper bound on the number of false positives in the 'x' most significant items
#' - TP: Lower bound on the number of true positives in the 'x' most significant items
#' - FDP: Upper bound on the proportion of false positives in the 'x' most significant items
#' - TP: Lower bound on the proportion of true positives in the 'x' most significant items
#'
#' Defaults to `c("TP", "FDP")`
#' @param all A logical value: should the bounds for all ordered subsets of `S` be returned? If `FALSE` (the default), only the bound for `S` is returned
#' @param ... Not used
#'
#' @return If `all` is `FALSE` (the default), only the value of the bound is returned. Otherwise, a `data.frame` is return, with |S| rows and 4 columns:
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
  bounds <- posthoc_bound2(p.values, S = S, thr = thr, lab = lab, what = what, all = all)
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
#' - FP: Upper bound on the number of false positives in the 'x' most significant items
#' - TP: Lower bound on the number of true positives in the 'x' most significant items
#' - FDP: Upper bound on the proportion of false positives in the 'x' most significant items
#' - TP: Lower bound on the proportion of true positives in the 'x' most significant items
#'
#' Defaults to `c("TP", "FDP")`
#' @param all A logical value: should the bounds for all ordered subsets of `S` be returned? If `FALSE` (the default), only the bound for `S` is returned
#'
#' @return If `all` is `FALSE` (the default), only the value of the bound is returned. Otherwise, a `data.frame` is return, with |S| rows and 4 columns:
#' - x: Number of most significant items selected
#' - label: Label for the procedure, typically of the form 'family(param)'
#' - bound: Value of the post hoc bound
#' - stat: Type of post hoc bound, as specified by argument `bound`.
#'
#' @export
#' @importFrom stats predict
#' @import sanssouci
posthoc_bound2 <- function(p.values, S = seq_along(p.values), thr = NULL, lab = NULL,
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
  table <- data.frame("Name" = c(), "# genes" = c(), "TP&ge;" = c(), "FDP&le;" = c(), check.names = FALSE)
  bioFun <- object$input$biologicalFunc
  if (class(bioFun)[1] == "list") {
    print("on passe ici")
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
      incProgress(1 / length(nameFunctions))
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
  return(table)
}

formatBounds <- function(max_FP, idxs = seq_len(max_FP), lab = NULL,
                         what = c("TP", "FDP"), all = FALSE) {
  stopifnot(length(max_FP) == length(idxs))
  max_FDP <- max_FP/idxs
  what0 <- c("FP", "TP", "FDP", "TDP")
  if (!all(what %in% what0)) {
    stop("Error in argument 'what': only the following statistics are supported: ", paste(what0, collapse = ", "))
  }
  if (length(max_FP) == 0) {
    if (all) {   # output should be empty (no subsets of positive size)
      mat <- matrix(NA, nrow = 0, ncol = 4)
      colnames(mat) <- c("x", "label", "stat", "bound")
      bounds <- as.data.frame(mat)
      return(bounds)
    } else {     # output should not be empty (number of FP in empty set is 0)
      idxs <- 0
      max_FP <- 0 # number of FP in empty set is 0
      max_FDP <- 0 # FDP in empty set is 0
    }
  }
  annot <- data.frame(x = idxs,
                      label = lab,
                      row.names = NULL)
  boundsList <- list(
    FP = cbind(annot, stat = "FP", bound = max_FP),
    TP = cbind(annot, stat = "TP", bound = idxs - max_FP),
    FDP = cbind(annot, stat = "FDP", bound = max_FDP),
    TDP = cbind(annot, stat = "TDP", bound = 1 - max_FDP))
  boundsList <- boundsList[what]
  if (!all) {
    boundsList <- lapply(boundsList, FUN = function(x) x[length(idxs), ])
  }
  bounds <- Reduce(rbind, boundsList)
  bounds
}
