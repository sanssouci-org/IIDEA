#' user interface of IIDEA
#'
#' @return UI structure of IIDEA
#' @export
#'
#' @import shiny
#' @importFrom plotly plotlyOutput
#' @import shinyBS
#' @importFrom shinyjs useShinyjs hidden
#' @importFrom DT DTOutput
#' @import htmlwidgets
app_ui <- function() {
  fluidPage(
    useShinyjs(),
    includeCSS(system.file("ShinyApps/www", "style.css", package = "IIDEA")),

    # Application title
    titlePanel(
      "IIDEA: Interactive Inference for Differential Expression Analyses"
    ),
    # Sidebar with panel
    sidebarLayout(
      sidebarPanel(
        tags$table(
          style = "width: 100%",
          tags$tr(

            tags$td(

              align = "center",
              htmlOutput("help")
            ),
            tags$td(
              align = "center",
              checkboxInput("checkboxDemo",
                            label = "Use public data",
                            value = TRUE
              )
            ),
            tags$td(
              align = "center",
              actionButton("buttonValidate", "Run!")
            )
          )
        ),

        selectInput("choiceTypeData",
                    label = "Choose the type of data",
                    choices = c(
                      "Microarrays" = "microarrays",
                      "bulk RNAseq" = "rnaseq"
                    )
        ),

        conditionalPanel(
          condition = "input.checkboxDemo",
          uiOutput("choiceGSEAUI")
        ),
        conditionalPanel(
          condition = "!input.checkboxDemo",
          fileInput("fileData",
                    label = p(
                      "Gene expression data matrix",
                      bsButton("QfileData",
                               label = "",
                               icon = icon("question"),
                               style = "info",
                               size = "extra-small"
                      ),

                      actionButton(
                        "resetInputData",
                        icon("trash")
                      )
                    ),
                    accept = ".csv"
          ),
          bsTooltip("QfileData",
                    "Upload a CSV file containing matrix with genes
                                 in rows and samples in column.
                                 Column names should be in (in {0, 1})",

                    "right",
                    options = list(container = "body"),
                    trigger = "hover"
          )
        ),
        conditionalPanel(
          condition = "!input.checkboxDemo",
          fileInput("fileLightData",
                    label = p(
                      "External p-values and log Fold Change matrix",
                      bsButton("QfileLightData",
                               label = "",
                               icon = icon("question"),
                               style = "info",
                               size = "extra-small"
                      ),

                      actionButton(
                        "resetInputLightData",
                        icon("trash")
                      )
                    ),
                    accept = ".csv"
          ),
          bsTooltip("QfileLightData",
                    "Upload a CSV file containing matrix with genes
                                 in rows and pvalues / Fold Change in columns.
                                 Column names should be 'fc', 'p.value'",

                    "right",
                    options = list(container = "body"),
                    trigger = "hover"
          )
        ),
        conditionalPanel(
          condition = "!input.checkboxDemo",

          # splitLayout(
          #   fileInput("fileAnnotation", label = p("Input gene annotation",
          #                                         bsButton("QfileAnnotation", label = "", icon = icon("question"), style = "info", size = "extra-small"))),

          fileInput("fileGroup",
                    label = p(
                      "Gene set matrix",
                      bsButton("QfileGroup",
                               label = "",
                               icon = icon("question"),
                               style = "info",
                               size = "extra-small"
                      ),

                      actionButton(
                        "resetInputGroup",
                        icon("trash")
                      )
                    ),
                    accept = ".csv"
          ),
          bsTooltip(
            id = "QfileGroup",
            title = "Upload a csv file containing matrix
                                 within nameGenes in line index. Binary vector
                                 composed this matrix for each gene set.",

            placement = "bottom",
            trigger = "hover",
            options = NULL
          ),

          downloadButton(
            "downloadExampleData",
            "Download example data set"
          )

        ),
        sliderInput("sliderConfLevel",
                    "Confidence level",
                    min = 0,
                    max = 100, value = 90, post = " %"

        ),
        checkboxInput("checkboxAdvancedParam",
                      label = p("Advanced parameters"),
                      # bsButton("Qparam",
                      #          label = "",
                      #          icon = icon("question"),
                      #          style = "info",
                      #          size = "extra-small")),
                      value = FALSE
        ),
        # bsTooltip(id = "Qparam",
        #           title = paste("Select parameters to implement permutation-based post hoc inference bounds for differential gene expression analysis, see dedicated ",
        #                           a("vignette.",
        #                             href = "https://sanssouci-org.github.io/sanssouci/articles/post-hoc_differential-expression.html")),
        #           trigger = c("click", "hover"),
        #           options = NULL),
        conditionalPanel(
          condition = "input.checkboxAdvancedParam",
          shinyjs::hidden(uiOutput("msgDegraded")),
          splitLayout(
            selectInput("alternative",
                        label = "Alternative",
                        choices = list(
                          "Two sided" = "two.sided",
                          "Less" = "less",
                          "Greater" = "greater"
                        ),
                        selected = "two.sided"
            ),
            numericInput("numB", label = "Number of permutations", value = 500, min = 10)
          ),
          # conditionalPanel(
          #   condition = "input.fileData",
          splitLayout(
            selectInput("refFamily",
                        label = "Reference family",
                        choices = list("Simes" = "Simes", "Beta" = "Beta"),
                        selected = "Simes"
            ),
            uiOutput("inputK") # )
            # )
          ),
          uiOutput("teststatUI")
        ),

        verbatimTextOutput("sorti"),
        conditionalPanel(
          condition = "input.buttonValidate != 0",
          tabsetPanel(
            id = "tabSelected",
            tabPanel("User selections",
                     value = 1,
                     uiOutput("OutQtableBounds"),
                     fluidRow(
                       column(
                         DTOutput("tableBounds"),
                         width = 12
                       )
                     ),

                     hidden(
                       downloadButton(
                         "downloadPHBTable",
                         "Download post hoc bound table"
                       )
                     )

            ),
            tabPanel("Gene sets",
                     value = 2,
                     uiOutput("OutQtableBoundsGroup"),
                     selectInput("buttonSEA",
                                 label = "Simultaneous Enrichment Analysis",
                                 choices = list(
                                   "All gene sets" = "nothing",
                                   "Significant for self-contained method" = "self",
                                   "Significant for competitive method" = "competitive"
                                 )
                     ),
                     uiOutput("errorMatch"),
                     DTOutput("tableBoundsGroup"),

                     downloadButton(
                       "downloadPHBTableGroup",
                       "Download post hoc bound table"
                     )

            )
          )
        ),
      ),

      # Main panel
      mainPanel(

        uiOutput("errorInput"),
        conditionalPanel(
          condition = "input.buttonValidate != 0",
          h2("Volcano plot",

             bsButton("Qparam1",
                      label = "",
                      icon = icon("question"), style = "info",
                      size = "extra-small"
             ),

             align = "center"
          ),
          bsPopover(
            id = "Qparam1",
            title = "VolcanoPlot",

            content = paste('Select genes by dragging
                                                  horizontal or vertical bars,
                                                  of using "box select" or
                                                  "lasso select" from the plot
                                                  menu. The table in the left
                                                  panel gives post-hoc bounds
                                                  for these selections.'),

            placement = "bottom",
            trigger = "hover",
            options = NULL
          ),
          flowLayout(
            selectInput("choiceYaxis",
                        label = "'y' axis label",
                        choices = list(
                          "p-values" = "pval",
                          "Adjusted p-values" = "adjPval",
                          "Number of false positves" = "thr"
                        ),
                        selected = "thr"
            ),
            checkboxInput("symetric",
                          label = "Symmetric fold change threshold",
                          value = FALSE
            ),
            uiOutput("msgURLds")
          )

        ),
        conditionalPanel(
          condition = "input.tabSelected==1",
          plotly::plotlyOutput("volcanoplotPosteriori",
                               height = "600px"
          ),
          fluidRow(
            shinyjs::hidden(actionButton(
              "resetCSV",
              "Reset Selections"
            )),
            shinyjs::hidden(downloadButton(
              "downloadData",
              "Download csv file with user selection"
            ))
          )
        ),
        conditionalPanel(
          condition = "input.tabSelected==2",
          uiOutput("errorBioMatrix"),
          plotly::plotlyOutput("volcanoplotPriori",
                               height = "600px"
          )

        )
      )
    ),
    p(em(
      "This interactive ",
      a("shiny", href = "https://shiny.rstudio.com", target = "_blank"),
      "application is developed by",

      a("Nicolas Enjalbert-Courrech",
        href = "https://nicolas-enjalbert.github.io/",
        target = "_blank"
      ),
      "and",
      a("Pierre Neuvial",
        href = "https://www.math.univ-toulouse.fr/~pneuvial/",
        target = "_blank"
      ),
      "for the R package ",
      a("sanssouci.",
        href = "https://sanssouci-org.github.io/sanssouci/",
        target = "_blank"
      ),
      "It implements permutation-based post hoc inference bounds for
       differential gene expression analysis, see dedicated ",
      a("vignette.",
        href = "https://sanssouci-org.github.io/sanssouci/articles/post-hoc_differential-expression.html",
        target = "_blank"
      ),
      "The ",
      a("source code",
        href = "https://github.com/sanssouci-org/sanssouci/tree/develop/inst/shiny-examples/volcano-plot",
        target = "_blank"
      ),
      "for this app is freely available. For any question, please file an",
      a("issue.",
        href = "https://github.com/sanssouci-org/sanssouci/issues",
        target = "_blank"
      )

    ))
  )
}
