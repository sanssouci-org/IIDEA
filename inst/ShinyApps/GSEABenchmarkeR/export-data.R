library(GSEABenchmarkeR)

# change work directory to create Example data set (work also with deployed
# version)
orwd <- getwd()
newwd <- dir(".", "GSEABenchmarkeR",
  all.files = TRUE, recursive = TRUE,
  full.names = TRUE, include.dirs = TRUE
)
setwd(paste0(getwd(), "/", newwd))

if (!dir.exists("express-data-set")) {
  dir.create("express-data-set")
}
if (!dir.exists("gene-set")) {
  dir.create("gene-set")
}

name_gse <- c(
  "GSE14762", "GSE15471", "GSE18842", "GSE19728", "GSE5281_EC",
  "GSE23878", "GSE7305", "GSE3467", "GSE9476", "GSE38666_epithelia"
)
data <- R.cache::memoizedCall(GSEABenchmarkeR::loadEData, "geo2kegg")
for (name in name_gse) {
  print(name)

  raw_data <- R.cache::memoizedCall(maPreproc, data[name])[[1]]

  saveRDS(raw_data, file = paste("express-data-set/", name, ".RDS", sep = ""))
  print("ok")
}


go_gs <- R.cache::memoizedCall(EnrichmentBrowser::getGenesets,
  org = "hsa", db = "go", onto = "BP",
  mode = "GO.db"
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
