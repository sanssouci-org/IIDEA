# IIDEA: Interactive Inference for Differential Expression Analyses

This package contains the shiny Application IIDEA. 

Run the online [IIDEA](https://shiny-iidea-sanssouci.apps.math.cnrs.fr/) shiny application.

## Install the package IIDEA 

``` r
remotes::install_github("sanssouci-org/IIDEA")
```

## Run offline IIDEA

To run offline IIDEA, please load the package IIDEA and then run the function `run_IIDEA` as 

``` r
IIDEA::run_IIDEA()
```
or 
``` r
library("IIDEA")
run_IIDEA()
```

To explore other example datasets, please perform the following functions

``` r 
load_microarray_datasets()
load_bulkRNAseq_datasets()
```

<!-- badges: start -->
[![R-CMD-check](https://github.com/sanssouci-org/IIDEA/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/sanssouci-org/IIDEA/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->
