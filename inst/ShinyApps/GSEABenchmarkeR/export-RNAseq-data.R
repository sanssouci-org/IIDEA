library(GSEABenchmarkeR)

# change work directory to create Example data set (work also with deployed version)
orwd <- getwd()
newwd <- dir(".", "GSEABenchmarkeR", all.files = TRUE, recursive=TRUE, full.names=TRUE, include.dirs=TRUE)
setwd(paste0(getwd(),'/', newwd))

if (!dir.exists("express-RNAseq-data-set")){
  dir.create("express-RNAseq-data-set")
}
if (!dir.exists("gene-set")){
  dir.create("gene-set")
}

nameGSE <- c("BRCA", "HNSC", "KICH", "KIRC", "KIRP",
             "LUSC", "PRAD", "STAD", "UCEC")
data <- R.cache::memoizedCall(GSEABenchmarkeR::loadEData,"tcga")
for (name in nameGSE){
  print(name)

  matrix <- SummarizedExperiment::assays(data[[name]])$exprs
  cats <- SummarizedExperiment::colData(data[[name]])
  ww <- match(cats$sample, base::colnames(matrix))
  categ <- cats$GROUP[ww]
  colnames(matrix) <- categ

  saveRDS(matrix, file=paste("express-RNAseq-data-set/",name,".RDS", sep = ""))
  print("ok")
}


go.gs <- R.cache::memoizedCall(EnrichmentBrowser::getGenesets,
                               org = "hsa", db = "go", onto = "BP", mode = "GO.db")

cleanGo.GS <- function(go.gs){
  for(i in names(go.gs)){
    if(length(go.gs[[i]])<10){
      go.gs[i] <- NULL
    }
  }
  return(go.gs)

}

go.gs <- R.cache::memoizedCall(cleanGo.GS, go.gs) # our func

saveRDS(go.gs, file="gene-set/go.gs.RDS")

print("go.gs done")

setwd(orwd)
