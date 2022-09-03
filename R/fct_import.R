#' Clean go.gs data sets. Not used
#'
#' @param go.gs
#'
#' @return a list of cleaned data set
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
  stop("Only ',', ';', '\\t' separators are allowed. Please adapte your csv file.")
}




#' Clean expression matrix gived by user
#'
#' @param matrixFunc a matrix for expression genes
#'
#' @return a list with the cleaned matrix, a text for error, a bool for error and the color of the error
#'
#' @details Clean the user matrix for IIDEA requirement. If the matrix does not
#' respect the requirement, an error is signaled. The color of the message means:
#' red: the matrix can be used
#' orange: the matrix doesnt respect but IIDEA update it and can be used
#' blue: everything is ok. Thanks ;)
#'
#' @export
cleanMatrix <- function(matrixFunc){
  text = ""
  boolValidation = TRUE
  color = "color:blue"

  for (i in 1:length(matrixFunc)){
    if(!is.numeric(matrixFunc[,i])){
      text <- "The matrix contains textual variable. We cannot use them."
      boolValidation <- FALSE
      color <- "color:red"
      return(list(data = matrixFunc, text = text, boolValidation = boolValidation, color = color))
    }
  }

  colnam <- colnames(matrixFunc) #Control colnames(matrix)
  valueColnam <- sort(unique(colnam))
  if(length(valueColnam)==2){
    if(any(valueColnam != c("0","1"))){
      colnam[ colnam == valueColnam[1] ] <- 0
      colnam[ colnam == valueColnam[2] ] <- 1
      colnames(matrixFunc) <- colnam
      text <- paste("The colnames of your matrix does not contains 0 or 1.We consider that", valueColnam[1] ,"becomes 0 and ", valueColnam[2]," becomes 1")
      boolValidation <- TRUE
      color <- "color:orange"
    }
  }else{
    boolValidation <- FALSE
    color <- "color:red"
    text <- "The column names of your data contains more (or less) than 2 categories. Please use {0, 1} for the colnames of your matrix."
  }
  return(list(data = matrixFunc, text = text, boolValidation = boolValidation, color = color))
}

#' Clean gene set matrix gived by user
#'
#' @param biofun a matrix for gene sets
#'
#' @return a list with the cleaned matrix, a text for error, a bool for error and the color of the error
#'
#' @details Clean the user matrix for IIDEA requirement. If the matrix does not
#' respect the requirement, an error is signaled. The color of the message means:
#' red: the matrix can be used
#' orange: the matrix doesnt respect but IIDEA update it and can be used
#' blue: everything is ok. Thanks ;)
#'
#' @export
cleanBiofun <- function(biofun){

  text = ""
  boolValidation = TRUE
  color = "color:blue"

  # if biofun cotains at least one textual variable
  for (i in 1:length(biofun)){
    if(!is.numeric(biofun[,i])){
      text <- "The gene set matrix contains textual variable. We cannot use them."
      boolValidation <- FALSE
      color <- "color:red"
      return(list(biofun = biofun, text = text, boolValidation = boolValidation, color = color))
    }
  }

  if(!all(biofun == 0 | biofun==1)){
    text <- "The gene set matrix contains other values than 0 or 1"
    boolValidation <- FALSE
    color <- "color:red"
  }

  return(list(biofun = biofun, text = text, boolValidation = boolValidation, color = color))
}

#' Check if genes in gene set match with a gene set list
#'
#' @param geneNames vector of genes names (vector of charactors)
#' @param biofun matrix for gene set
#'
#' @return
#' @export
#'
#' #' @details Check if there is at least one gene in gene set match with genes
#' in expression matrix.  The color of the message means:
#' green: non genes match
#' blue: everything is ok. Thanks ;)
#'
#' @examples
matchMatrixBiofun <- function(geneNames, biofun){
  text = ""
  boolValidation = TRUE
  color = "color:blue"

  mm <- match(rownames(biofun), geneNames)
  biofun <- biofun[mm, ]
  if(all(is.na(mm))){
    text = "None of the lines of the gene set matrix correspond to the lines of the gene expression data matrix."
    boolValidation = FALSE
    color = "color:green"
  }

  return(list(biofun = biofun, text = text, boolValidation = boolValidation, color = color))

}

