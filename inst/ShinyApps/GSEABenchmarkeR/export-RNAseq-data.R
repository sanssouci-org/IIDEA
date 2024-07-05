library(GSEABenchmarkeR)

# change work directory to create Example data set (work also with deployed
# version)
orwd <- getwd()
newwd <- dir(".", "GSEABenchmarkeR",
  all.files = TRUE, recursive = TRUE,
  full.names = TRUE, include.dirs = TRUE
)
setwd(paste0(getwd(), "/", newwd))

if (!dir.exists("express-RNAseq-data-set")) {
  dir.create("express-RNAseq-data-set")
}
if (!dir.exists("gene-set")) {
  dir.create("gene-set")
}

name_gse <- c(
  "BRCA", "HNSC", "KICH", "KIRC", "KIRP",
  "LUSC", "PRAD", "STAD", "UCEC"
)
data <- R.cache::memoizedCall(GSEABenchmarkeR::loadEData, "tcga")
for (name in name_gse) {
  print(name)

  matrix <- SummarizedExperiment::assays(data[[name]])$exprs
  cats <- SummarizedExperiment::colData(data[[name]])
  ww <- match(cats$sample, base::colnames(matrix))
  categ <- cats$GROUP[ww]
  colnames(matrix) <- categ

  saveRDS(matrix, file = paste("express-RNAseq-data-set/", name, ".RDS",
                               sep = ""))
  print("ok")
}


go_gs <- R.cache::memoizedCall(EnrichmentBrowser::getGenesets,
  org = "hsa", db = "go", onto = "BP", mode = "GO.db"
)

clean_go_gs <- function(go_gs) {
  for (i in names(go_gs)) {
    if (length(go_gs[[i]]) < 10) {
      go_gs[i] <- NULL
    }
  }
  return(go_gs)
}

go_gs <- R.cache::memoizedCall(clean_go_gs, go_gs) # our func

saveRDS(go_gs, file = "gene-set/go.gs.RDS")

print("go.gs done")

setwd(orwd)
