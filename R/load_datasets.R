#' Install GSEAbenchmarkeR dataset from the GEO2KEGG microarray compendium to be used in IIDEA
#'
#' @param names a vector of character of GSEAbenchmarkeR expression data
#' compendium.  Should be a subset of c("GSE14762", "GSE15471", "GSE18842",
#' "GSE19728", "GSE5281_EC","GSE23878", "GSE7305", "GSE3467","GSE9476",
#' "GSE38666_epithelia")
#'
#' @return NULL
#'
#' @importFrom GSEABenchmarkeR loadEData maPreproc
#' @importFrom R.cache memoizedCall
#' @importFrom EnrichmentBrowser getGenesets
#' @importFrom SummarizedExperiment assays colData
load_microarray_datasets <- function(names = c(
  "GSE14762", "GSE15471",
  "GSE18842", "GSE19728",
  "GSE5281_EC",
  "GSE23878", "GSE7305", "GSE3467",
  "GSE9476", "GSE38666_epithelia"
)) {

  names_possible <- c(
    "GSE14762", "GSE15471",
    "GSE18842", "GSE19728",
    "GSE5281_EC",
    "GSE23878", "GSE7305", "GSE3467",
    "GSE9476", "GSE38666_epithelia"
  )

  if (!all(names %in% names_possible)) {
    stop(sprintf(
      "The argument 'names' should be a subset of c('%s')",
      paste0(names_possible, collapse = "', '")
    ))
  }

  path <- system.file("ShinyApps", package = "IIDEA")
  if (!dir.exists(file.path(path, "GSEABenchmarkeR"))) {
    dir.create(file.path(path, "GSEABenchmarkeR"))
  }
  path <- file.path(path, "GSEABenchmarkeR")
  path_express_data <- "express-data-set"
  if (!dir.exists(file.path(path, path_express_data))) {
    dir.create(file.path(path, path_express_data))
  }
  path_gene_set <- "gene-set"
  if (!dir.exists(file.path(path, path_gene_set))) {
    dir.create(file.path(path, path_gene_set))
  }


  pattern <- "(.*) \\((.*)\\).RDS"
  filenames <- list.files(file.path(path, path_express_data), pattern = pattern)

  exist_rds <- gsub(pattern, "\\2", filenames)
  names <- setdiff(names, exist_rds)

  data <- R.cache::memoizedCall(GSEABenchmarkeR::loadEData, "geo2kegg")
  for (name in names) {
    print(name)
    rawData <- R.cache::memoizedCall(GSEABenchmarkeR::maPreproc, data[name])[[1]]
    matrix <- SummarizedExperiment::assays(rawData)$exprs

    cats <- SummarizedExperiment::colData(rawData)
    ww <- match(cats$Sample, base::colnames(matrix))
    categ <- cats$GROUP[ww]
    colnames(matrix) <- categ

    lID<- paste(rawData@metadata$experimentData@other$disease, " (", (name), ")", sep = "")

    saveRDS(matrix,
      file = file.path(
        path,
        path_express_data,
        paste(lID, ".RDS", sep = "")
      )
    )
  }

  if (!file.exists(file.path(path, path_gene_set, "go.gs.RDS"))) {
    cleanGo.GS <- function(go.gs) {
      for (i in names(go.gs)) {
        if (length(go.gs[[i]]) < 10) {
          go.gs[i] <- NULL
        }
      }
      return(go.gs)
    }

    go.gs <- R.cache::memoizedCall(EnrichmentBrowser::getGenesets,
      org = "hsa", db = "go", onto = "BP",
      mode = "GO.db"
    )

    go.gs <- R.cache::memoizedCall(cleanGo.GS, go.gs) # our func

    saveRDS(go.gs, file = file.path(path, path_gene_set, "go.gs.RDS"))
  }
}

#' Install GSEAbenchmarkeR dataset from the TCGA RNA-seq compendium to be used in IIDEA
#'
#' @param names vector of character of GSEAbenchmarkeR expression data
#' compendium.  Should be a subset of c("BRCA", "HNSC", "KICH", "KIRC", "KIRP",
#' "LUSC", "PRAD", "STAD", "UCEC")
#'
#' @return NULL
#'
#' @importFrom GSEABenchmarkeR loadEData
#' @importFrom SummarizedExperiment assays colData
#' @importFrom R.cache memoizedCall
#' @importFrom EnrichmentBrowser getGenesets
load_bulkRNAseq_datasets <- function(names = c(
  "BRCA", "HNSC", "KICH", "KIRC", "KIRP",
  "LUSC", "PRAD", "STAD", "UCEC"
)) {
  names_possible <- c(
    "BRCA", "HNSC", "KICH", "KIRC", "KIRP",
    "LUSC", "PRAD", "STAD", "UCEC"
  )

  if (!all(names %in% names_possible)) {
    stop(sprintf(
      "The argument 'names' should be a subset of c('%s')",
      paste0(names_possible, collapse = "', '")
    ))
  }


  path <- system.file("ShinyApps", package = "IIDEA")
  if (!dir.exists(file.path(path, "GSEABenchmarkeR"))) {
    dir.create(file.path(path, "GSEABenchmarkeR"))
  }
  path <- file.path(path, "GSEABenchmarkeR")
  path_express_data <- "express-RNAseq-data-set"
  if (!dir.exists(file.path(path, path_express_data))) {
    dir.create(file.path(path, path_express_data))
  }
  path_gene_set <- "gene-set"
  if (!dir.exists(file.path(path, path_gene_set))) {
    dir.create(file.path(path, path_gene_set))
  }

  # pattern <- "(.*)_\\((.*)\\).RDS"
  pattern <- "(.*).RDS"
  filenames <- list.files(file.path(path, path_express_data), pattern = pattern)

  exist_rds <- gsub(pattern, "\\2", filenames)
  names <- setdiff(names, exist_rds)

  data <- R.cache::memoizedCall(GSEABenchmarkeR::loadEData, "tcga")
  for (name in names) {
    matrix <- SummarizedExperiment::assays(data[[name]])$exprs
    cats <- SummarizedExperiment::colData(data[[name]])
    ww <- match(cats$sample, base::colnames(matrix))
    categ <- cats$GROUP[ww]
    colnames(matrix) <- categ

    CPM <- matrix / colSums(matrix) * 1e6
    row_maxs <- matrixStats::rowMaxs(CPM)
    ww <- which(row_maxs < 10)
    row_quantiles <- matrixStats::rowQuantiles(log(1 + CPM),
                                               prob = 0.75
    )
    ww <- which(row_quantiles < log(1 + 5))
    matrix <- log(1 + CPM[-ww, ])

    saveRDS(matrix,
      file = file.path(
        path,
        path_express_data,
        paste(name, ".RDS", sep = "")
      )
    )
  }


  cleanGo.GS <- function(go.gs) {
    for (i in names(go.gs)) {
      if (length(go.gs[[i]]) < 10) {
        go.gs[i] <- NULL
      }
    }
    return(go.gs)
  }

  go.gs <- R.cache::memoizedCall(EnrichmentBrowser::getGenesets,
    org = "hsa", db = "go", onto = "BP", mode = "GO.db"
  )

  go.gs <- R.cache::memoizedCall(cleanGo.GS, go.gs) # our func

  saveRDS(go.gs, file = file.path(path, path_gene_set, "go.gs.RDS"))

  print("go.gs done")
}
