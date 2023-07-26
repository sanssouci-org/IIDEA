#' Clean go.gs data sets. Not used
#'
#' @param go.gs data frame from go.gs
#'
#' @return a list of cleaned data set
#' @export
#' @importFrom shiny incProgress
cleanGo.GS <- function(go.gs) {
  for (i in names(go.gs)) {
    incProgress(1 / length(go.gs))
    if (length(go.gs[[i]]) < 6) {
      go.gs[i] <- NULL
    }
  }
  return(go.gs)
}

#' Read multiple csv configuration
#'
#' @param file a character of the csv file path
#' @param ... other parameter for read.csv function
#'
#' @return a dataframe containing the dataset from the csv file
#' @importFrom utils read.csv
#' @export
readCSV_sep <- function(file, ...) {
  df <- read.csv(file, sep = c(","), ...)
  if (dim(df)[2] > 1) {
    return(df)
  }
  df <- read.csv(file, sep = c(";"), ...)
  if (dim(df)[2] > 1) {
    return(df)
  }
  df <- read.csv(file, sep = c("\t"), ...)
  if (dim(df)[2] > 1) {
    return(df)
  }
  stop("Only ',', ';', '\\t' separators are allowed.
       Please adapte your csv file.")
}




#' Clean expression matrix gived by user
#'
#' @param matrixFunc a matrix for expression genes
#'
#' @return a list with the cleaned matrix, a text for error,
#' a bool for error and the color of the error
#'
#' @details Clean the user matrix for IIDEA requirement. If the matrix does not
#' respect the requirement, an error is signaled.
#' The color of the message means:
#' red: the matrix can be used
#' orange: the matrix doesnt respect but IIDEA update it and can be used
#' blue: everything is ok. Thanks ;)
#'
#' @export
cleanMatrix <- function(matrixFunc) {
  text <- ""
  boolValidation <- TRUE
  color <- "color:blue"

  for (i in seq_along(matrixFunc)) {
    if (!is.numeric(matrixFunc[, i])) {
      text <- "The matrix contains textual variable. We cannot use them."
      boolValidation <- FALSE
      color <- "color:red"
      return(list(
        data = matrixFunc, text = text,
        boolValidation = boolValidation, color = color
      ))
    }
  }

  colnam <- colnames(matrixFunc) # Control colnames(matrix)
  valueColnam <- sort(unique(colnam))
  if (length(valueColnam) == 2) {
    if (any(valueColnam != c("0", "1"))) {
      colnam[colnam == valueColnam[1]] <- 0
      colnam[colnam == valueColnam[2]] <- 1
      colnames(matrixFunc) <- colnam
      text <- paste(
        "The colnames of your matrix does not contains 0 or 1.
                    We consider that", valueColnam[1],
        "becomes 0 and ", valueColnam[2], " becomes 1"
      )
      boolValidation <- TRUE
      color <- "color:orange"
    }
  } else {
    boolValidation <- FALSE
    color <- "color:red"
    text <- "The column names of your data contains more (or less) than
    2 categories. Please use {0, 1} for the colnames of your matrix."
  }
  return(list(
    data = matrixFunc, text = text, boolValidation = boolValidation,
    color = color
  ))
}

#' Clean gene set matrix gived by user
#'
#' @param biofun a matrix for gene sets
#'
#' @return a list with the cleaned matrix, a text for error,
#' a bool for error and the color of the error
#'
#' @details Clean the user matrix for IIDEA requirement. If the matrix does not
#' respect the requirement, an error is signaled.
#'  The color of the message means:
#' red: the matrix can be used
#' orange: the matrix doesnt respect but IIDEA update it and can be used
#' blue: everything is ok. Thanks ;)
#'
#' @export
cleanBiofun <- function(biofun) {
  text <- ""
  boolValidation <- TRUE
  color <- "color:blue"

  # if biofun cotains at least one textual variable
  for (i in seq_along(biofun)) {
    if (!is.numeric(biofun[, i])) {
      text <- "The gene set matrix contains textual variable.
      We cannot use them."
      boolValidation <- FALSE
      color <- "color:red"
      return(list(
        biofun = biofun, text = text, boolValidation = boolValidation,
        color = color
      ))
    }
  }

  if (!all(biofun == 0 | biofun == 1)) {
    text <- "The gene set matrix contains other values than 0 or 1"
    boolValidation <- FALSE
    color <- "color:red"
  }

  return(list(
    biofun = biofun, text = text, boolValidation = boolValidation,
    color = color
  ))
}

#' Check if genes in gene set match with a gene set list
#'
#' @param geneNames vector of genes names (vector of charactors)
#' @param biofun matrix for gene set
#'
#' @return list for warnings
#'
#' @details Check if there is at least one gene in gene set match with genes
#' in expression matrix.  The color of the message means:
#' green: non genes match
#' blue: everything is ok.
#'
#' @export
matchMatrixBiofun <- function(geneNames, biofun) {
  text <- ""
  boolValidation <- TRUE
  color <- "color:blue"

  mm <- match(rownames(biofun), geneNames)
  biofun <- biofun[mm, ]
  if (all(is.na(mm))) {
    text <- "None of the lines of the gene set matrix correspond to the lines of
    the gene expression data matrix."
    boolValidation <- FALSE
    color <- "color:green"
  }

  return(list(
    biofun = biofun, text = text, boolValidation = boolValidation,
    color = color
  ))
}


#' Build gene set for BLCA RNAseq data set from sanssouci.data
#'
#' @param RNAseq_blca matrix
#'
#' @return matrix
#' @export
#' @import org.Hs.eg.db
RNAseq_blca_GO <- function(RNAseq_blca){
  # BiocManager::install("org.Hs.eg.db")
  library("org.Hs.eg.db")

  # data set-specific code to perform ad hoc gene selection
  # (to avoid having too many GO terms)
  dat <- RNAseq_blca
  categ <- ifelse(colnames(dat) == "III", 1, 0) # map to 0/1
  dex <- data.frame(rowWilcoxonTests(dat, categ))
  probeNames <- rownames(dex)
  nbProbes <- length(probeNames)

  pval <- dex[["p.value"]]
  adjp <- p.adjust(pval, method = "BH")
  #selected <- probeNames[which(adjp < 0.05)]
  selected <- probeNames[which(pval < 0.0001)]
  length(selected)
  # end data set-specific code

  keytypes(org.Hs.eg.db)

  cols <- c("SYMBOL", "GENENAME")
  sel <- select(org.Hs.eg.db, keys=selected, columns=cols, keytype="SYMBOL")
  head(sel)
  dim(sel)

  sel <- select(org.Hs.eg.db, keys=selected, columns="GO", keytype="SYMBOL")
  head(sel)
  dim(sel)
  sel <- subset(sel, EVIDENCE == "TAS")
  dim(sel)

  GOs <- unique(sel$GO)
  rev_sel <- select(org.Hs.eg.db, keys=GOs, columns=cols, keytype="GO")
  head(rev_sel)
  rev_sel <- subset(rev_sel, EVIDENCE == "TAS")

  head(rev_sel)
  dim(rev_sel)
  table(rev_sel$ONTOLOGY)

  # map to matrix with genes x go_terms
  gnames <- rownames(dat)
  tap <- tapply(rev_sel$SYMBOL, INDEX = rev_sel$GO, FUN = function(x){
    y <- numeric(length(gnames))
    y[match(x, gnames)] <- 1
    y
  })
  go <- Reduce(cbind, tap)
  rownames(go) <- gnames
  colnames(go) <- names(tap)
  dim(go)

  expr_ALL_GO <- go[, colSums(go) >= 3]

  return(expr_ALL_GO)
}

