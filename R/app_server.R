#' Server part of IIDEA
#'
#' @param input automatic inputs of the shiny server
#' @param output automatic output
#' @param session automatic session
#'
#' @return a server instance
#' @export
#' @import sanssouci.data
#' @import sanssouci
#' @import shiny
#' @import shinyjs
#' @importFrom matrixStats rowMaxs rowQuantiles
#' @importFrom plotly event_data plotlyProxyInvoke plot_ly add_markers event_register
#' config toWebGL plotlyProxy layout
#' @importFrom stats p.adjust
#' @importFrom utils zip str write.csv
#' @importFrom stringr str_remove_all
#' @importFrom SummarizedExperiment assays colData
#' @importFrom dplyr select filter
app_server <- function(input, output, session) {

  data(expr_ALL, package = "sanssouci.data", envir = environment())
  data(expr_ALL_GO, package = "sanssouci.data", envir = environment())
  data(RNAseq_blca, package = "sanssouci.data", envir = environment())
  data(RNAseq_blca_GO, package = "sanssouci.data", envir = environment())

  # size of input data sets
  options(shiny.maxRequestSize = 1024^3)
  output$help <- renderUI({
    a("IIDEA help page",
      href = "https://sanssouci-org.github.io/sanssouci/articles/IIDEA.html",
      target = "_blank"
    )

  })

  ###################
  # Loading data
  ###################

  ## button run

  isolate({

    disable("buttonValidate")
  }) # while geo2kegg is not loaded, user cannot "run" #Initialization

  observeEvent(object_I(), { # object is loaded, user can "run"
    enable("buttonValidate")

  })

  ## input for example data sets

  namesExampleFile <- reactive({

    if(req(input$choiceTypeData) == "microarrays"){

      filenames <- (list.files("inst/ShinyApps/GSEABenchmarkeR/express-data-set",
                               pattern = "*.RDS", full.names = TRUE
      ))

      if (length(filenames) == 0) {
        return(NULL)
      }
      ldf <- lapply(filenames[2:length(filenames)], readRDS)
      lID <- sapply(ldf, function(l) {
        l@metadata$dataId
      })
      names(lID) <- paste(sapply(ldf, function(l) {
        l@metadata$experimentData@other$disease
      }), " (", (lID), ")", sep = "")
      return(lID)}
    else if (req(input$choiceTypeData) == "rnaseq"){
      filenames <- (list.files("inst/ShinyApps/GSEABenchmarkeR/express-RNAseq-data-set",
                               pattern = "*.RDS", full.names = TRUE
      ))

      if (length(filenames) == 0) {
        return(NULL)
      }

      pattern <- "inst/ShinyApps/GSEABenchmarkeR/express-RNAseq-data-set/(.*).RDS"
      lID <- sapply(filenames[-1], function(filename) {gsub(pattern, "\\1", filename)})
      names(lID) = lID
      return(lID)

    }
  }) # get names of data sets


  output$choiceGSEAUI <- renderUI({ # create input

    req(input$choiceTypeData)
    if (req(input$choiceTypeData) == "microarrays") {
      selectInput("choiceGSEA",
                  label = "Choose a gene data set",
                  choices = c(
                    "Leukemia (ALL): BCR/ABL mutated vs wild type" = "OurData",
                    namesExampleFile()
                  )
      )
    } else if (req(input$choiceTypeData) == "rnaseq") {
      selectInput("choiceGSEA",
                  label = "Choose a gene data set",
                  choices = c("Urothelial Bladder Carcinoma: stage II vs stage III" = "BLCA",
                              namesExampleFile()
                  )
      )
    }

  })

  ## Download our example data set, loaded from sanssouci.data

  ### load data sets in the button : save in three csv file containing
  # in a zip file
  output$downloadExampleData <- downloadHandler(
    filename = function() {
      paste("ExampleData", "zip", sep = ".")
    },
    content = function(fname) {
      exampleData <- exampleData(req(input$choiceTypeData))
      fs <- c()
      tmpdir <- tempdir()
      setwd(tempdir())
      paths <- c("expressData.csv", "biologicalFunction.csv", "volcanoData.csv")
      write.csv(exampleData$matrix, paths[1])
      write.csv(exampleData$biologicalFunc, paths[2])
      write.csv(exampleData$degrade, paths[3])
      zip(zipfile = fname, files = paths)
    },
    contentType = "application/zip"
  )

  ## Loading user inputs

  ### Loading expression gene matrix
  fileData <- reactiveVal(NULL) # Initialisation
  observe({ # when a new csv file is in input
    newValue <- req(input$fileData)
    fileData(newValue)
  })

  observeEvent(input$resetInputData, {
    # to delete input file (clicking bin incon)

    fileData(NULL) # delete serveur variable
    data(NULL)
    reset("fileData") # delete UI variable (in fileInput)
  })

  ### Loading expression gene matrix
  fileLightData <- reactiveVal(NULL) # Initialisation
  observe({ # when a new csv file is in input
    newValue <- req(input$fileLightData)
    fileLightData(newValue)
  })

  observeEvent(input$resetInputLightData, {
    # to delete input file (clicking bin incon)

    fileLightData(NULL) # delete serveur variable
    data(NULL)
    reset("fileLightData") # delete UI variable (in fileInput)
  })

  ## Loading gene set matrix : same structure as before
  fileGroup <- reactiveVal(NULL)
  observe({
    newValue <- input$fileGroup
    fileGroup(newValue)
  })
  observeEvent(input$resetInputGroup, {
    fileGroup(NULL)
    data(NULL)
    reset("fileGroup")
  })

  pvaluesVP <- reactiveVal(NULL)
  logFCVP <- reactiveVal(NULL)

  ## Cleaning data
  object_I <-
    reactive({
      # eventReactive(input$buttonValidate, {
      withProgress(value = 0, message = "Upload Data... ", {
        if (input$checkboxDemo) { # if example data set

          if (req(input$choiceTypeData) == "microarrays") {
            if (req(input$choiceGSEA) == "OurData" | req(input$choiceGSEA) == "BLCA") {
              # cleaning for data from sanssouci.data
              # print("Etape 1")
              setProgress(value = 0.4, detail = "sanssouci data set ...")
              matrix <- expr_ALL # read data from sanssouci.data
              # print("Etape 2")
              ### cleaning categories
              cat <- colnames(matrix)
              categ <- rep(1, length(cat))
              categ[which(cat == "NEG")] <- 0
              # print("Etape 3")

              object <- SansSouci(Y = as.matrix(matrix), groups = as.numeric(categ)) # create SansSouci object
              object$input$geneNames <- rownames(matrix)
              # print("Etape 4")
              setProgress(value = 0.7, detail = "Preparation of gene set data ...  ")

              bioFun <- expr_ALL_GO
              stopifnot(nrow(bioFun) == nrow(matrix)) ## sanity check: dimensions
              mm <- match(base::rownames(bioFun), base::rownames(matrix))
              stopifnot(!any(is.na(mm)))
              object$input$biologicalFunc <- bioFun[mm, ]
              rm(bioFun)
              # print("Etape 5")
              object$bool$validation <- TRUE
              object$bool$degrade <- FALSE
              rm(matrix)
              rm(categ)
            } else {
              # cleaning data set from GSEA data set
              setProgress(value = 0.4, detail = "GSEA data set ...")
              rawData <- readRDS(paste("inst/ShinyApps/GSEABenchmarkeR/express-data-set/", input$choiceGSEA, ".RDS", sep = ""))

              matrix <- assays(rawData)$exprs

              cats <- colData(rawData)
              ww <- match(cats$Sample, base::colnames(matrix))
              categ <- cats$GROUP[ww]
              object <- SansSouci(Y = matrix, groups = as.numeric(categ))
              setProgress(value = 0.7, detail = "GSEA data set ...")

              object$input$geneNames <- base::rownames(matrix)

              object$input$biologicalFunc <- readRDS("inst/ShinyApps/GSEABenchmarkeR/gene-set/go.gs.RDS")
              # On a laisse sous forme de liste car on a adapte les fonctions qui en ont besoin. Plus rapide qu'en la transformant en matrice binaire
              object$bool$url <- rawData@metadata$experimentData@url

              object$bool$validation <- TRUE
              object$bool$degrade <- FALSE
            }
          } else if (req(input$choiceTypeData) == "rnaseq") {
            if (req(input$choiceGSEA) == "BLCA" | req(input$choiceGSEA) == "OurData") {
              # cleaning for data from sanssouci.data
              setProgress(value = 0.4, detail = "sanssouci data set ...")
              matrix <- RNAseq_blca # read data from sanssouci.data

              CPM <- matrix / colSums(matrix) * 1e6
              # plot(density(rowMaxs(log(1 + CPM))))
              row_quantiles <- matrixStats::rowQuantiles(log(1 + CPM), prob = 0.75)
              ww <- which(row_quantiles < log(1 + 5))
              matrix <- log(1 + CPM[-ww, ])

              ### cleaning categories
              cat <- colnames(matrix)
              categ <- rep(1, length(cat))
              categ[which(cat == "II")] <- 0

              object <- SansSouci(Y = as.matrix(matrix), groups = as.numeric(categ)) # create SansSouci object
              object$input$geneNames <- rownames(matrix)

              setProgress(value = 0.7, detail = "Preparation of gene set data ...  ")

              bioFun <- RNAseq_blca_GO # read data from sanssouci.data
              ## sanity check: all genes are annotated
              mm <- match(base::rownames(matrix), base::rownames(bioFun))
              stopifnot(!any(is.na(mm)))
              bioFun <- bioFun[mm, ]
              stopifnot(identical(base::rownames(matrix), base::rownames(bioFun)))

              # make sure at least 3 genes for each GO *after filtering*
              bioFun <- bioFun[, colSums(bioFun) >= 3]
              str(bioFun)

              object$input$biologicalFunc <- bioFun
              rm(bioFun)
              object$bool$validation <- TRUE
              object$bool$degrade <- FALSE
              rm(matrix)
              rm(categ)
            } else {
              setProgress(value = 0.4, detail = "GSEA data set ...")
              rawData <- readRDS(paste("inst/ShinyApps/GSEABenchmarkeR/express-RNAseq-data-set/", input$choiceGSEA, ".RDS", sep = ""))
              # rawData <- R.cache::memoizedCall(maPreproc,geo2kegg()[input$choiceGSEA])[[1]]

              matrix <- rawData

              CPM <- matrix / colSums(matrix) * 1e6
              # plot(density(rowMaxs(log(1 + CPM))))
              row_maxs <- matrixStats::rowMaxs(CPM)
              ww <- which(row_maxs < 10)
              row_quantiles <- matrixStats::rowQuantiles(log(1 + CPM), prob = 0.75)
              ww <- which(row_quantiles < log(1 + 5))
              matrix <- log(1 + CPM[-ww, ])

              categ <- colnames(matrix)
              object <- SansSouci(Y = matrix, groups = as.numeric(categ))
              setProgress(value = 0.7, detail = "GSEA data set ...")

              object$input$geneNames <- base::rownames(matrix)

              object$input$biologicalFunc <- readRDS("inst/ShinyApps/GSEABenchmarkeR/gene-set/go.gs.RDS")
              # object$bool$url <- rawData@metadata$experimentData@url

              object$bool$validation <- TRUE
              object$bool$degrade <- FALSE

            }
          }
        } else {
          # if user use his own data
          # req(input$fileData, input$fileLightData)
          # req(fileData(), fileLightData())
          req(isTruthy(fileData()) || isTruthy(fileLightData()))
          # file <- req(fileData())

          is.expressionMatrix <- !is.null((fileData()))
          is.VPMatrix <- !is.null((fileLightData()))
          # print("Etape 5")
          setProgress(value = 0.1, detail = "Read csv ...")

          if(is.expressionMatrix){
            matrix <- readCSV_sep(file = fileData()$datapath, row.names = 1, check.names = FALSE)
            # print("Etape 6")
            setProgress(value = 0.4, detail = "Clean data set ...")
            clean <- cleanMatrix(matrix) # see cleanMatrix function :
            ### give specific issue for non available matrix
            # print("Etape 7")
            ### if matrix is not available, matrix == NULL allows to block the calibration JER and the printing
            if (clean$boolValidation) { # if matrix is ok
              object <- SansSouci(Y = as.matrix(clean$data), groups = as.numeric(colnames(clean$data)))

              object$input$geneNames <- base::rownames(clean$data)
            } else { # if isn't
              input <- NULL
              parameters <- NULL
              output <- NULL

              object <- structure(
                list(
                  input = input,
                  parameters = parameters,
                  output = output
                ),
                class = "SansSouci"
              )
            }

            object$bool$validation <- clean$boolValidation
            object$bool$matrix.color <- clean$color # color of error message

            object$bool$matrix.text <- clean$text # error message
          } else {
            setProgress(value = 0.2, detail = "Test data ...")

            input <- NULL
            parameters <- NULL
            output <- NULL

            object <- structure(
              list(
                input = input,
                parameters = parameters,
                output = output
              ),
              class = "SansSouci"
            )
          }

          if (is.VPMatrix){
            VPmatrix <- readCSV_sep(file = fileLightData()$datapath,
                                    row.names = 1, check.names = FALSE)
            setProgress(value = 0.2, detail = "Test data ...")
            pvaluesVP((VPmatrix[["p.value"]]))
            logFCVP(VPmatrix[["fc"]])
            setProgress(value = 0.4, detail = "Add pvalues and logFC ...")
            if(!is.expressionMatrix){
              m <- dim(VPmatrix)[1]
              input <- list(m = m, geneNames = base::rownames(VPmatrix))
              alpha <- 0.1
              parameters <- list(
                B = 0,
                family = "Simes", k = m
              )

              output <- list(
                p.value = VPmatrix[["p.value"]],
                estimate = VPmatrix[["fc"]]
              )

              object$input = input
              object$parameters = parameters
              object$output = output
              object$bool$validation <- TRUE
            }

          }

          object$bool$degrade <- !is.expressionMatrix & is.VPMatrix


          setProgress(value = 0.7, detail = "Preparation of gene set data ...")

          ## cleaning of gene set matrix
          fileGroup <- fileGroup()
          if (!is.null(fileGroup)) {
            setProgress(value = 0.75, detail = "Read gene set data ...")
            T1 <- Sys.time()

            bioFun <- readCSV_sep(
              file = fileGroup$datapath, row.names = 1,
              check.names = FALSE
            )
            T2 <- Sys.time()

            setProgress(value = 0.8, detail = "Cleaning of gene set data ...")
            cleanBio <- cleanBiofun(bioFun) # cleaning and message error

            rm(bioFun)

            object$bool$bioFun.color <- cleanBio$color
            object$bool$bioFun.text <- cleanBio$text


            setProgress(value = 0.9, detail = "Matching ...")


            matchBio <- matchMatrixBiofun(
              geneNames = object$input$geneNames,
              biofun = cleanBio$biofun
            )
            # verification compatibility between the two matrices
            if (matchBio$boolValidation & cleanBio$boolValidation) {

              object$input$biologicalFunc <- as.matrix(matchBio$biofun)
            } else { # if not ok
              object$input$biologicalFunc <- NULL
            }
            object$bool$match.color <- matchBio$color
            object$bool$match.text <- matchBio$text

            rm(matchBio)
            rm(cleanBio)
          }

          rm(matrix)
        }
        setProgress(value = 1, detail = "Done")
        # print("Etape exit object_I()")
        return(object)
      })
    })

  data <- reactiveVal()
  observe({
    data(req(object_I()))
    # print("update data")
  })




  urlDataSet <- eventReactive(input$buttonValidate, {
    # give link to description of geo2kegg data sets

    req(input$choiceGSEA)
    req(object_I()$bool$url)
    return(a("URL of data set description",
             href = object_I()$bool$url,
             target = "_blank"
    ))

  })
  output$msgURLds <- renderUI({
    req(urlDataSet)
    tagList(urlDataSet())
  })


  ## different error messages for non compliant matrix
  # express gene matrix error
  output$errorInput <- renderUI({
    tags$span(
      style = req(object_I()$bool$matrix.color),
      paste(req(object_I()$bool$matrix.text))
    )
  })
  # gene set matrix error
  output$errorBioMatrix <- renderUI({
    tags$span(
      style = req(object_I()$bool$bioFun.color),
      paste(req(object_I()$bool$bioFun.text))
    )
  })
  # no common gene b/w expression gene matrix and gene set matrix
  output$errorMatch <- renderUI({
    tags$span(
      style = req(object_I()$bool$match.color),
      paste(req(object_I()$bool$match.text))
    )
  })


  output$watch <- renderPrint({
    tableCSV()
  })



  ###################
  # Parameters
  ###################

  # Confidance alpha
  alpha <- reactiveVal(0.1) # Initialization

  observeEvent(input$buttonValidate, {
    # When Run is clicked : we get the input value

    newValue <- req(1 - input$sliderConfLevel / 100)
    alpha(newValue)
  })

  # number of permutation
  numB <- reactiveVal(500) # Initialization

  observeEvent(input$buttonValidate, {
    # When Run is clicked : we get the input value

    newValue <- req(input$numB)
    numB(newValue)
  })

  # Test statistic
  output$teststatUI <- renderUI({ # create input
    # req(input$choiceTypeData)

    selectInput("teststat",
                label = "Choose a test statistic for the calibration",
                choices = switch(req(input$choiceTypeData),
                                 microarrays = c(
                                   "Welch test" = rowWelchTests
                                 ),
                                 rnaseq = c(
                                   "Wilcoxon test" = rowWilcoxonTests
                                 )
                )
    )
  })


  rowTestFUN <- reactiveVal()
  observeEvent(object_I(), {
    req(input$choiceTypeData)
    newValue <- switch(input$choiceTypeData,
                       microarrays = rowWelchTests,
                       rnaseq = rowWilcoxonTests
    )
  })

  # Attention ici le test n'est pas mis a jour
  observeEvent(input$teststat, {
    req(input$choiceTypeData)
    newValue <- req(input$teststat)
    rowTestFUN(newValue)
  })

  # reference family
  refFamily <- reactiveVal("Simes") # Initialization

  observeEvent(input$buttonValidate, {
    # When Run is clicked : we get the input value

    newValue <- req(input$refFamily)
    refFamily(newValue)
  })

  # alternative hypothesis

  alternative <- reactiveVal("two.sided") # Initialization
  observeEvent(input$buttonValidate, {
    # When Run is clicked : we get the input value

    newValue <- req(input$alternative)
    alternative(newValue)
  })

  # parameter K
  ## dynamic input (need matrixChosen())
  output$inputK <- renderUI({
    req(object_I())
    # print("entree output$inputtK")
    numericInput("valueK",
                 label = "K (size of reference family)",
                 value = numKI(),
                 min = 1,
                 max = nrow(req(object_I()$input$Y))
    )
    # print("sortie output$inputtK")
  })

  ## numKI() is used to intiate the printed input valueK
  numKI <- reactiveVal()
  observe({ # Initialization, if refFamily == 'Beta' or not
    req(object_I())
    # print("entree numKI")
    newValue <- ifelse(input$refFamily == "Beta",
                       round(2 * req(nrow(object_I()$input$Y)) / 100),
                       req(nrow(object_I()$input$Y))
    )
    numKI(newValue)
    # print("sorti numKI")
  })

  ## numK is the parameters choosen by users and use in server side
  numK <- reactiveVal()

  # if the parameters are not activated, input$valueK doesn't exist
  # and therefore we don't have the right result...
  # why ? CalibrateJER should handle if K is null, right?
  observeEvent(object_I(), {
    # print("entree numK")
    req(object_I()) # when object_I() is change, INITIALISATION of numK()
    newValue <- req(nrow(object_I()$input$Y))
    numK(newValue)
    # print("sortie numK")
  })
  # When Run is clicked : we get the input value.
  # Here an issue : if, advanced parameters is not opened, input$valueK == NULL
  observeEvent(input$buttonValidate, {
    newValue <- req(input$valueK)
    numK(newValue)
  })



  # If light matrix is available

  output$msgLight <- renderUI({
    tags$span(
      style = "color:grey",
      paste(
        "A matrix containing p-values and fold change is detected.",
        "Thus, you cannot change the following advanced parameters:\n",
        "Reference family = 'Simes',\n K = ", length(object_I()$output$p.value)
      )
    )

  })

  observe({
    req(object_I())
    if (object_I()$bool$degrade) {

      show("msgLight")


      hide("alternative")
      hide("numB")
      hide("refFamily")
      hide("inputK")
    } else {
      hide("msgLight")

      show("alternative")
      show("numB")
      show("refFamily")
      show("inputK")
    }
  })


  ###################
  # Calibration
  ###################

  ## JER calibration on available expression gene matrix
  # observe({
  observeEvent(input$buttonValidate, {
    req(data()$bool$validation)

    is.expressionMatrix <- !is.null((fileData()))
    is.VPMatrix <- !is.null((fileLightData()))

    if (is.expressionMatrix | input$checkboxDemo) { # non light version

      withProgress(value = 0, message = "Perform calibration ... ", {
        incProgress(amount = 0.3)
        rTF <- req(rowTestFUN())
        if (!is.function(class(rTF))) {
          rTF <- req(eval(parse(text = rowTestFUN()), envir = environment(rowWelchTests)))
        }
        t1 <- Sys.time()
        object <- fit(data(),
                      alpha = req(alpha()),
                      B = numB(),
                      rowTestFUN = rTF,
                      alternative = alternative(),
                      family = refFamily(),
                      K = numK()
        )

        t2 <- Sys.time()

        if (!is.VPMatrix | input$checkboxDemo){
          pvaluesVP((pValues(object)))
          logFCVP(foldChanges(object))
        }
        setProgress(value = 0.7, detail = "Done")
      })

    } else { # light version
      # matrices with pval and fc are in the object since its creation

      m <- nHyp(data())
      thr <- t_linear(alpha(), seq_len(m), m)
      # force using of Simes and k=m # IMPORT FROM FUNCTION.R

      object <- data()
      object$output$thr <- thr
      object$output$lambda <- alpha()
      object$parameters$alpha <- alpha()
    }
    # calcul des logp et adjp
    object$output$logp <- -log10(pValues(object))
    object$output$adjp <- p.adjust(pValues(object), method = "BH")
    data(object)
  })

  logpvaluesVP <- reactiveVal(NULL)
  pvaluesAdjVP <- reactiveVal(NULL)
  observe({
    req(pvaluesVP())
    # print("entree pvaluesupdate")
    logpvaluesVP(-log10(pvaluesVP()))
    newValue <- p.adjust(pvaluesVP(), method = "BH")
    pvaluesAdjVP(newValue)
    # print("sortie pvaluesupdate")
  })


  ###################
  # Threshold
  ###################


  # vertical contains the displacement value of all movable objects
  # (in our case just the thresholds)
  vertical <- reactive({
    event_data("plotly_relayout", source = "A")
  })


  # Updating user movement
  ## threshold logfc right
  xint <- reactiveVal(0.5)

  observeEvent(vertical()[["shapes[0].x0"]], {
    # activate when the object 0 change

    if (input$symetric) {
      newValue <- vertical()[["shapes[0].x0"]]
      xint(abs(newValue))
      xint2(-abs(newValue))
    } else {
      newValue <- vertical()[["shapes[0].x0"]]
      if(newValue <xint2()){
        exchange <- xint2()
        xint2(newValue)
        newValue <- exchange
      }
      xint(newValue)
    }
  })

  ## threshold logfc left
  xint2 <- reactiveVal(-0.5)

  observeEvent(vertical()[["shapes[2].x0"]], {
    # activate when the object 2 change

    if (input$symetric) {
      newValue <- vertical()[["shapes[2].x0"]]
      xint(abs(newValue))
      xint2(-abs(newValue))
    } else {
      newValue <- vertical()[["shapes[2].x0"]]
      if(newValue > xint()){
        exchange <- xint()
        xint(newValue)
        newValue <- exchange
      }
      xint2(newValue)
    }
  })

  ## threshold pval
  yint <- reactiveVal()
  observeEvent(pvaluesVP(), {
    # initialization : adjusted p-value == 0.1
    req(data())
    # print("entree yint()")
    min0.1 <- which.min(abs(pvaluesAdjVP() - 0.1))
    y0.1 <- logpvaluesVP()[[min0.1]]
    yint(y0.1)
    # print("sortie yint()")
  })
  observeEvent(vertical()[["shapes[1].y0"]], {
    # when user change threshold on plotly
    p <- 10^(-vertical()[["shapes[1].y0"]])

    y_sel <- which((pvaluesVP() <= p)) # selected by  p-value
    newValue <- Inf
    if (length(y_sel) > 0) {
      newValue <- min(logpvaluesVP()[y_sel])
      # threshold on the log(p-value) scale

    }
    yint(newValue)
  })

  ###################
  # Selecting gene with thresholds
  ###################

  # selected genes for server calcuation
  selectedGenes <- reactive({
    req(logFCVP(),logpvaluesVP())
    # data.frame(x = req(logFCVP()), y = req(pvaluesVP()))
    ## gene selections
    # print("entree selectedGenes")
    sel1 <- which(logpvaluesVP() >= yint() &
                    logFCVP() >= xint()) # upper right
    sel2 <- which(logpvaluesVP() >= yint() &
                    logFCVP() <= xint2()) # upper left

    sel12 <- sort(union(sel1, sel2)) # both
    # print("sortie selectedGenes")
    return(list(sel1 = sel1, sel2 = sel2, sel12 = sel12))

  })



  # matrix p-val & fc of selected genes by thresholds
  # => to plot red points on volcano plot
  selected_points <- reactive({
    # print("selected_points")
    # print(length(selectedGenes()$sel12))
    list(
      x = logFCVP()[selectedGenes()$sel12],
      y = logpvaluesVP()[selectedGenes()$sel12]
    )
  })

  ###################
  # Selecting gene by lasso or box
  ###################

  # reactive variable containing value of lasso/box selecting
  d <- reactive({
    event_data("plotly_selected", source = "A")
  })

  # list of selected genes (list of name rows)
  manuelSelected <- reactive({
    req(d())
    d()[["pointNumber"]][which(d()[["curveNumber"]] == 0)] + 1
  })


  ###################
  # Calculate post hob bounds
  ###################

  # POST Hoc bound (PHB) for thresholds selection
  TP_FDP <- reactive({
    req(selectedGenes())

    ## post hoc bounds in selections
    req(data())
    req(thresholds(data()))
    # print("entree TP_FDP")
    n12 <- length(selectedGenes()$sel12)
    pred <- predict(
      object = data(), S = selectedGenes()$sel12,
      what = c("TP", "FDP")
    )
    # print("sortie TP_FDP")
    return(list(
      n12 = n12, TP12 = pred["TP"], FDP12 = pred["FDP"]
    ))
  })

  # calculate PHB for lasso box selection
  calcBoundSelection <- reactive({ #
    req(manuelSelected())
    req(thresholds(data()))
    # print("calcBoundSelection")
    c(n = length(manuelSelected()), predict(data(),
                                            S = manuelSelected(),
                                            what = c("TP", "FDP")
    ))

  })



  ###################
  # Post Hoc bound table reactive variable
  ###################

  ## Initialization of variable of post hoc bound table (PHB table)
  tableResult <- reactiveVal(data.frame(
    Selection = c("Threshold selection")
  )) # Initialization

  # PHB table for only selected genes by thresholds
  # reactive : TP_FDP() <- selectedGenes() <- {xint(), xint2(), yint()} <-
  # vertical() <- user threshold moving
  baseTable <- reactive({
    req(TP_FDP())
    data.frame(
      `Selection` = c("Threshold selection"),
      "# genes" = c(TP_FDP()$n12),
      "TP\u2265" = as.integer(c(TP_FDP()$TP12)),
      "FDP\u2264" = c(round(TP_FDP()$FDP12, 2)),
      check.names = FALSE, row.names = NULL
    )
  })

  # updating PHB table when thresholds change
  observeEvent(TP_FDP(), { # When threshold change

    print("entree tableresults")
    print(tableResult())
    bottomTable <- tableResult() %>% # keep box/lasso selection
      filter(`Selection` != "Threshold selection")
    # remove row named "Threshold selection"
    upperTable <- baseTable()
    # take new value of PHB for the new threshold selection

    newValue <- rbind(upperTable, bottomTable)
    tableResult(newValue)
    # print("sortie tableresults")
  })

  # Add PHB for lasso and box selection
  ## d() is the reactive variable for box and lasso selection
  observeEvent(d(), { # When user selects a new group of points
    req(calcBoundSelection())
    # print("entree d()")

    vectorGene <- names(pValues(data())[manuelSelected()])
    # list of gene contained in gene selection
    url <- UrlStringdbGrah(vectorGene)
    # construction of link to StringDB interaction graph
    n <- dim(tableResult())[1]
    newValue <- rbind(tableResult(), c(
      paste('<a target="_blank" href="', url, '" >User selection ',
            n, "</a>", sep = ""),

      calcBoundSelection()["n"],
      calcBoundSelection()["TP"],
      round(calcBoundSelection()["FDP"], 2)
    ))
    tableResult(newValue)
    # print("sortie d()")
  })

  # To clean gene selection from lasso/box selection
  ## keep only threshold selection
  observeEvent(input$resetCSV, { # to clean printed table
    newValue <- baseTable()
    tableResult(newValue)
  })

  # if data changes, PHB table is cleaned
  observeEvent(data(), { # to clean printed table
    newValue <- baseTable()
    tableResult(newValue)
  })


  ###################
  # Post Hoc bound table outputs
  ###################

  # reactive popify to explain PHB table
  output$OutQtableBounds <- renderUI({
    tab <- tableResult()
    msg <- "The table below prints post hoc bounds for user selections. "
    if (nrow(tab) > 0) {
      msg <- paste(
        msg, "For example, the selection called",
        tab[1, 1],
        "contains at least",
        tab[1, "TP\u2265"],

        " true positives (TP) and its False Discovery Proportion (FDP)
        is less than",

        tab[1, "FDP\u2264"]
      )
    }
    popify(
      el = bsButton("QtableBounds",
                    label = "", icon = icon("question"),
                    style = "info", size = "extra-small"
      ),

      title = "Data",
      content = msg,
      trigger = "hover"
    )
  })

  # output for PHB table
  output$tableBounds <- renderDT(
    {
      req(TP_FDP())
      tableResult()
    },

    selection = list(mode = "single", selectable = -(1))
    # can select only one row and not the first one

    ,
    escape = FALSE # to print url link to open stringDB graph
  )

  ###################
  # Prepare objects for plotly volcano plots
  ###################

  # value for adjusted p-values yaxis
  lineAdjp <- reactive({
    req(pvaluesAdjVP())
    listLog <- c()
    # print("entree linAdjp")
    for (i in c(0.5, 0.25, 0.1, 0.05, 0.025, 0.01, 0.001, 0.0001)) {
      # selected line on yaxis /!\ if you change it you have to change
      # in yaxis()
      min05 <- which.min(abs(pvaluesAdjVP() - i))
      # to take the gene with the nearest adjpvalue to i
      y05 <- logpvaluesVP()[[min05]]
      # take the equivalent in logpvalue to print it on VP

      listLog <- c(listLog, y05)
    }
    # print("sortie linAdjp")
    return(listLog)
  })

  # value for 'NUmber of false discoverie' yaxis : optimize y line
  thr_yaxis <- reactive({
    req(alpha()) # if parameters change (not only alpha)
    req(data())
    req(thresholds(data()))
    req(logpvaluesVP())
    tya <- thrYaxis(thr = thresholds(data()), maxlogp = max(logpvaluesVP()))
    return(tya)
  })

  # reactive variable containing values for yaxis depending users choice

  yaxis <- reactive({
    req(data())
    f <- list(
      size = 14,
      color = "#000000"
    )
    yaxis <- switch(input$choiceYaxis,
                    "pval" = list(
                      title = "p-value (-log[10] scale)",
                      titlefont = f
                    ),
                    "adjPval" = list(
                      title = "Adjusted p-value (-log[10] scale)",
                      titlefont = f,
                      autotick = FALSE,
                      tickmode = "array",
                      tickvals = lineAdjp(),

                      ticktext = c(0.5, 0.25, 0.1, 0.05, 0.025, 0.01, 0.001, 0.0001)
                      # be careful of changing of lingAdjp()

                    ),
                    "thr" = list(
                      title = "Maximal number of false positives",
                      titlefont = f,
                      autotick = FALSE,
                      tickmode = "array",
                      tickvals = req(thr_yaxis()$pvalue),
                      ticktext = req(thr_yaxis()$num)
                    )
    )
    return(yaxis)
  })

  # reactive values for threshold (used for selecting genes)
  thrLine <- reactive({
    req(data())
    list(
      list( # right logFC threshold vertical()[["shapes[0].x0"]]
        type = "line",
        line = list(color = "orange", dash = "dot"),
        x0 = xint(),
        x1 = xint(),
        y0 = 0,
        y1 = 1,
        yref = "paper"
      ),
      list( # p-val threshold vertical()[["shapes[1].y0"]]
        type = "line",
        line = list(color = "orange", dash = "dot"),
        x0 = 0,
        x1 = 1,
        y0 = yint(),
        y1 = yint(),
        xref = "paper"
      ),
      list( # left logFC threshold vertical()[["shapes[2].x0"]]
        type = "line",
        line = list(color = "orange", dash = "dot"),
        x0 = xint2(),
        x1 = xint2(),
        y0 = 0,
        y1 = 1,
        yref = "paper"
      )
    )
  })

  ###################
  # OUTPUT object of plotly volcano plot VP1
  ###################

  # Intialisation of volcano plot for threshold selection (part 1)
  posteriori <- reactive({
    req(data()$bool$validation)
    req(data())

    req(logFCVP(), logpvaluesVP())
    setProgress(value = 0.9, detail = "posteriori Reactive ... ")
    f <- list(
      size = 14,
      color = "#000000"
    )
    lte <- "\u2264"
    gte <- "\u2265"
    p <- plot_ly(data.frame(x = req(logFCVP()), y = req(logpvaluesVP())),
                 x = ~x, y = ~y,
                 marker = list(
                   size = 2,
                   color = "grey"
                 ),
                 name = "unselected",
                 type = "scattergl", mode = "markers", # make points
                 source = "A" # name plotly::plot

    ) %>%
      add_markers(
        x = req(selected_points()$x), y = req(selected_points()$y),
        # add at the begening selected points in red

        marker = list(
          color = "red",
          size = 6
        ),
        name = "selected"
      ) %>%
      layout(
        xaxis = list(title = "Fold change (log scale)", titlefont = f),

        yaxis = isolate(yaxis()),
        # isolate is used not to reactive this reactive variable
        # (improve global reactivity : see bellow for the changing of yaxis)

        title = "",
        shapes = isolate(thrLine()), # same as before
        dragmode = "select", # initialise lasso/box selection by default
        showlegend = FALSE
      ) %>%

      event_register("plotly_selecting") %>% # to save point selected by
      # lasso:box selection

      config(editable = TRUE) %>%
      toWebGL() # to go faster on web
    return(p)
  })



  output$volcanoplotPosteriori <- renderPlotly({
    withProgress(message = "Posterio plot ...", value = 0, {
      p <- posteriori()
      shiny::setProgress(value = 1, detail = "Done")
      return(p)
    })
  })



  ###################
  # Change of VP1
  ###################

  # This part aims of changing only hight layers of the VP1 plotly object.
  # Indeed, posteriori() is only reactive if df() change.
  # Each 'observeEvent' represents an user change.
  # The VP1 is updated with the function plorlyProxy (choose the plotly OUTPUT)
  # and the plotlyProxyInvoke() function to define the action.


  # when threshold moved and red selected points change
  observeEvent(vertical(), {
    plotlyProxy("volcanoplotPosteriori", session) %>%
      plotlyProxyInvoke("deleteTraces", 1) %>% # delete the previous red points
      plotlyProxyInvoke(
        "addTraces", # add red points from selected_points()
        list(
          x = unname(selected_points()$x),
          y = unname(selected_points()$y),
          type = "scatter",
          mode = "markers",
          line = list(color = "red"),
          name = "selected",
          showlegend = TRUE
        ), 1 # the rank on the stack
      )
  })

  # # Yaxis change
  observeEvent(
    {
      input$choiceYaxis # user selection
      yaxis()
    } # update yaxis values when parameters change (alpha, ...)
    ,
    { # when we choose a different y axis
      plotlyProxy("volcanoplotPosteriori", session) %>%

        plotlyProxyInvoke("relayout", list(yaxis = yaxis()))
      # relayout to update values

      # a stack is not used in layout, only one value is possible to yaxis
    }
  )

  # When user changes orange threshold to select red points
  observeEvent(vertical(), { # for moving of orange threshold
    plotlyProxy("volcanoplotPosteriori", session) %>%
      plotlyProxyInvoke(
        "relayout",
        list(shapes = thrLine())
      ) # same as before
  })

  # on the stack TRACES of VP1, the layer 0 is all the point (unselected),
  # the layer 1 is the red points
  # the layer 2 is the blue points



  # when user want to watch again a user selection by box or lasso
  # tableBounds_rows_selected : the user selection on the output tableBounds
  ## take the name of the user selection
  userDTselectPost <- reactive({
    href <- tableResult()[input$tableBounds_rows_selected, "Selection"]
    str_remove_all(str_remove_all(href, "<a(.*?)>"), "(</a>)")
  })

  # take the list of gene selected by user
  selectionUserRe <- reactive({
    vect <- tableCSV()[, req(userDTselectPost())]
    sel <- which(vect == 1)
    list(sel = sel)
  })

  # update plotly VP1
  observeEvent(userDTselectPost(), {

    if (length(userDTselectPost()) == 1) {
      # if a row is selected (length == 1 because can select only one row)
      plotlyProxy("volcanoplotPosteriori", session) %>%
        plotlyProxyInvoke("deleteTraces", 2)
      # first we delete previous blue points :
      # if there are not bleu points, none points are deleted

      plotlyProxy("volcanoplotPosteriori", session) %>%
        # second, we trace new blue points

        plotlyProxyInvoke(
          "addTraces",
          list(
            x = unname(logFCVP()[selectionUserRe()$sel]),
            y = unname(logpvaluesVP()[selectionUserRe()$sel]),
            type = "scattergl",
            mode = "markers",
            line = list(color = "blue"),
            name = userDTselectPost()
          ), 2 # put it in layer 2
        )
    } else { # If none, two or more row are selected
      plotlyProxy("volcanoplotPosteriori", session) %>%
        plotlyProxyInvoke("deleteTraces", 2)
    }
  })


  ###################
  # Download list of selected gene from box/lasso
  ###################

  tableCSV <- reactiveVal() # creation of reactive variable

  # initialize table when data() change
  observeEvent(
    {
      input$resetCSV
      data()$input$Y
      input$buttonValidate
    },
    {
      req(data()$input$Y)
      req(selectedGenes())
      req(data()$input$geneNames)
      vecteur <- rep(0, dim(data()$input$Y)[1])
      vecteur[selectedGenes()$sel12] <- 1

      req(length(vecteur) == length(data()$input$geneNames))
      print(length(vecteur))
      print(length(data()$input$geneNames))

      newValue <- data.frame(
        Thresholds_selection = vecteur,
        row.names = req(data()$input$geneNames)
      )
      # a dataframe with rowname without features

      tableCSV(newValue)
    }
  )

  # When user want to change thresholds
  observeEvent(selectedGenes(), {
    req(data()$input$Y)
    req(selectedGenes())
    req(tableCSV())
    vecteur <- rep(0, dim(data()$input$Y)[1])
    vecteur[selectedGenes()$sel12] <- 1

    # print("avant 1229")
    rigthTable <- tableCSV() %>%
      dplyr::select(-Thresholds_selection)
    # print("apres 1232")
    df <- cbind(

      data.frame(
        Thresholds_selection = vecteur,
        row.names = req(data()$input$geneNames)
      ),

      rigthTable
    )
    tableCSV(df)
  })

  # When user selected ne gene set from lasso/box [d() activate]
  observeEvent(d(), {
    req(data()$input$Y)
    req(tableCSV())
    vecteur <- rep(0, dim(data()$input$Y)[1])
    vecteur[manuelSelected()] <- 1
    nameCol <- colnames(tableCSV())
    df <- cbind(tableCSV(), selection2 = vecteur)
    colnames(df) <- c(nameCol, paste("User selection", length(nameCol)))
    tableCSV(df)
  })



  # download binary matrix containing list of gene in user selections
  output$downloadData <- downloadHandler( # download csv of user selection
    filename = function() {
      tag <- format(Sys.time(), "%Y-%M-%d_%H-%m-%S")
      sprintf("volcano-plot_gene-selections_%s.csv", tag)
    },
    content = function(file) {
      write.csv(tableCSV(), file)
    }
  )

  # download PHB table
  output$downloadPHBTable <- downloadHandler(
    filename = function() {
      tag <- format(Sys.time(), "%Y-%M-%d_%H-%m-%S")
      sprintf("volcano-plot_bounds_%s.csv", tag)
    },
    content = function(file) {
      table <- tableResult()
      table$`Selection` <- str_remove_all(str_remove_all(
        table$`Selection`,
        "<a(.*?)>"
      ), "(</a>)")
      write.csv(table, file)
    }
  )

  # while buttonValidate is not clicked, these buttons are hidden
  observeEvent(input$buttonValidate, {
    show("downloadPHBTable")
    show("resetCSV")
    show("downloadData")
  })




  ###################
  # PART 2 : gene set analyses
  ###################

  # PHB for all gene sets
  tableBoundsGroup <- reactive({
    withProgress(message = "tableBoundsGroup", {
      T1 <- Sys.time()
      req(thresholds(req(data())))
      req(data()$input$biologicalFunc)

      table <- boundGroup2(req(data()))
      T2 <- Sys.time()
    })
    return(table)
  })


  # calculate vounds for all features to compare with each gene sets
  boundsW <- reactive({ # calculate bounds for all features

    req(pValues(data()))
    req(thresholds(data()))
    c(n = length(nHyp(data())), predict(object = data()))
  })


  # dowload csv file containing PHB table of gene sets

  output$downloadPHBTableGroup <- downloadHandler(
    # download csv of user selection

    filename = function() {
      tag <- format(Sys.time(), "%Y-%M-%d_%H-%m-%S")
      sprintf("gene-set_bounds_%s.csv", tag)
    },
    content = function(file) {
      table <- tableBoundsGroup()
      table$Name <- str_remove_all(str_remove_all(
        table$Name,
        "<a(.*?)>"
      ), "(</a>)")
      write.csv(table, file)
    }
  )

  # If user choose SEA alternative
  filteredTableBoundsGroup <- reactive({
    req(input$buttonSEA)
    req(tableBoundsGroup())
    if (input$buttonSEA == "competitive") {
      table <- tableBoundsGroup()
      sel <- which(table[["FDP\u2264;"]] < boundsW()["FDP"])
      newValue <- table[sel, ]
      return(newValue)
    } else if (input$buttonSEA == "self") {
      table <- tableBoundsGroup()
      sel <- which(table[["TP\u2265;"]] > 0)
      return(table[sel, ])
    } else {
      return(tableBoundsGroup())
    }
  })

  # reactive popify to explain PHB table
  output$OutQtableBoundsGroup <- renderUI({
    req(tableBoundsGroup())
    popify(

      el = bsButton("QtableBoundsGroup",
                    label = "",
                    icon = icon("question"), style = "info",
                    size = "extra-small"
      ),
      title = "Data",
      content = paste(

        "This table prints your post-hoc bounds for your gene sets.",
        "For example, the selection called",
        tableBoundsGroup()[1, "Name"],
        "contains at leat",
        tableBoundsGroup()[1, "TP\u2265"],

        " true positives (TP) and its False Discovery Proportion (FDP)
        is less than ",

        round(tableBoundsGroup()[1, "FDP\u2264"] * 100, 2), "%"
      ),
      trigger = "hover"
    )
  })


  output$tableBoundsGroup <- renderDT(
    {
      table <- filteredTableBoundsGroup()
      table[["FDP\u2264;"]] <- round(table[["FDP\u2264;"]], 2)

      table
    },
    selection = "single",
    escape = FALSE,
    options = list(scrollX = TRUE)
  )

  # name of gene set selected by user
  userDTselectPrio <- reactive({
    req(filteredTableBoundsGroup())

    href <- filteredTableBoundsGroup()[
      input$tableBoundsGroup_rows_selected,
      "Name"
    ]

    name <- str_remove_all(str_remove_all(href, "<a(.*?)>"), "(</a>)")
    return(name)
  })

  # list of genes selected by user
  selectionGroup <- reactive({
    req(data()$input$biologicalFunc)

    group <- req(userDTselectPrio())
    bioFun <- data()$input$biologicalFunc

    if (inherits(bioFun, "list")) {

      ids <- bioFun[[group]]
    } else {
      ids <- which(bioFun[, group] == 1)
    }
    list(sel = ids)
  })


  # 'reactive" plot : as the plot before, this one should be activate once.
  # See posteriori() for details

  # VP2
  priori <- reactive({
    req(data())
    req(logFCVP())
    req(logpvaluesVP())
    f <- list(
      size = 14,
      color = "#000000"
    )
    lte <- "\u2264"
    gte <- "\u2265"

    plot_ly(data.frame(x = logFCVP(), y = logpvaluesVP()),
            x = ~x, y = ~y,
            marker = list(
              size = 2,
              color = "grey"
            ),
            name = "genes",
            type = "scattergl", mode = "markers", source = "B",
            text = data()$input$geneNames,

            customdata = paste0(
              "http://www.ensembl.org/Homo_sapiens/Gene/Summary?g=",
              data()$input$geneNames
            )

    ) %>%
      layout(
        showlegend = TRUE,
        xaxis = list(title = "Fold change (log scale)", titlefont = f),
        # yaxis = isolate(yaxis()),
        title = "",
        dragmode = "select"
      ) %>%
      onRender("
                  function(el) {
                      el.on('plotly_click', function(d) {
                          var url = d.points[0].customdata;
                          window.open(url);
                      });
                  }
              ") %>%
      event_register("plotly_selecting") %>%
      config(editable = TRUE) %>%
      toWebGL()
  })

  # # output of priori()
  output$volcanoplotPriori <- renderPlotly({
    withProgress(message = "Plot", {
      p <- priori()
      shiny::setProgress(value = 1, detail = "Done")
      return(p)
    })
  })


  # when yaxis changes

  observeEvent(
    {
      input$choiceYaxis
      yaxis()
    },
    { # when we choose a different y axis
      plotlyProxy("volcanoplotPriori", session) %>%
        plotlyProxyInvoke("relayout", list(yaxis = yaxis()))
    }
  )


  # when user select a gen set to print it on VP2

  # here, there are no red points : stack is composed of 0 :points ;
  # 1 : blue points (gene set selected)

  observeEvent(userDTselectPrio(), {
    if (length(userDTselectPrio()) == 1) {
      plotlyProxy("volcanoplotPriori", session) %>%
        plotlyProxyInvoke("deleteTraces", 1)
      plotlyProxy("volcanoplotPriori", session) %>%
        plotlyProxyInvoke(
          "addTraces",
          list(
            x = unname(logFCVP()[selectionGroup()$sel]),
            y = unname(logpvaluesVP()[selectionGroup()$sel]),
            type = "scattergl",
            mode = "markers",
            line = list(color = "blue"),
            name = userDTselectPrio()
          )
        )
    } else {
      plotlyProxy("volcanoplotPriori", session) %>%
        plotlyProxyInvoke("deleteTraces", 1)
    }
  })

}
