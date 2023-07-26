utils::globalVariables(c("expr_ALL", "expr_ALL_GO"))
utils::globalVariables(c("RNAseq_blca", "RNAseq_blca_GO"))

#' Create example data set
#'
#' @param type character data type
#'
#' @details Create example data set from the package sanssouci.data.
#' This dataframe aims at exporting in IIDEA (not to be used by IIDEA as public
#' data set)
#'
#' @export
#'
#' @import sanssouci.data
#' @import sanssouci
#'
#' @examples
#' exampleData(type = "microarrays")
exampleData <- function(type) {

  if (type ==  "microarrays"){
    if (!exists("expr_ALL")) {
      data(expr_ALL, package = "sanssouci.data", envir = environment())
    }
    if (!exists("expr_ALL_GO")) {
      data(expr_ALL_GO, package = "sanssouci.data", envir = environment())
    }
    #### expression matrix
    data <- list()
    data$matrix <- expr_ALL

    categ <- colnames(data$matrix)
    data$categ <- rep(1, length(categ))
    data$categ[which(categ == "NEG")] <- 0

    colnames(data$matrix) <- data$categ

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
  } else if (type == "rnaseq") {
    if (!exists("RNAseq_blca")) {
      data(RNAseq_blca, package = "sanssouci.data", envir = environment())
    }
    if (!exists("RNAseq_blca_GO")) {
      data(RNAseq_blca_GO, package = "sanssouci.data", envir = environment())
    }
    
    #### expression matrix
    data <- list()
    data$matrix <- RNAseq_blca

    categ <- colnames(data$matrix)
    data$categ <- rep(1, length(categ))
    data$categ[which(categ == "II")] <- 0

    colnames(data$matrix) <- data$categ

    #### degraded matrix
    dex <- rowWilcoxonTests(data$matrix, data$categ)
    data$degrade <- data.frame("p.value" = dex[["p.value"]], "fc" = dex$estimate)
    rm(dex)

    # gene set matrix
    bioFun <- RNAseq_blca_GO
    mm <- match(base::rownames(bioFun), base::rownames(data$matrix))
    data$biologicalFunc <- bioFun[mm, ]
    rm(bioFun)
    return(data)
  }
}


#' Create an URL for string-db gene network from a vector of gene names.
#'
#' @param vector a vector of gene names use for the gene networks
#' from string-db.
#'
#' @return a character, the url for string-db.org
#'
#' @export
UrlStringdbGrah <- function(vector) {
  vector[2:length(vector)] <- paste0("%0d", vector[2:length(vector)])
  url <- paste("https://string-db.org/api/image/network?identifiers=",
               paste(vector, collapse = ""), "&species=9606", sep = "")
  return(url)
}

#' Create a html balise for Gene Ontology description on ebi.ac.uk
#'
#' @param name charactor of Gene Ontology name.
#'
#' @return charactor of html balise
#' @importFrom stringr str_extract_all
#' @export
#'
addUrlLink <- function(name) {
  if (grepl("GO:\\d+", name)) {
    url <- paste("https://www.ebi.ac.uk/QuickGO/term/",
                 str_extract_all(name, "GO:\\d+"), sep = "")
    url <- paste('<a target="_blanck" href="', url, '" >',
                 name, "</a>", sep = "")
    return(url)
  } else {
    return(name)
  }
}

#' Create a new y-axis for Volcano plot based on calibration threshold.
#'
#' @param thr vector for calibration threshold
#' @param maxlogp numeric the maximum of -log(p-values)
#'
#' @return data.frame with the new labels and breaks for the Volcano Plot y-axis
#' @export
thrYaxis <- function(thr, maxlogp) {
  if (maxlogp == Inf) {
    maxlogp <- 16
  }
  df1 <- data.frame(num = seq_along(thr) - 1, pvalue = -log10(thr))
  df2 <- data.frame(df1[c(1), ])
  valeurTest <- df2[c(dim(df2)[1]), "pvalue"]
  for (i in seq_len(dim(df1)[1])) {
    mod <- if (df1[i, "num"] < 100) {
      1
    } else if (df1[i, "num"] < 500) {
      10
    } else if (df1[i, "num"] < 1000) {
      50
    } else {
      100
    }
    if ((valeurTest - df1[i, "pvalue"] > 0.3 * maxlogp / 12.5)
        && (df1[i, "num"] %% mod == 0)) {
      df2 <- rbind(df2, (df1[i, ]))
      valeurTest <- df1[i, "pvalue"]
    }
  }
  return(df2)
}
