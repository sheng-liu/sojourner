
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Sojourner: Statistical Analysis of Single Molecule Trajectories

<!-- badges: start -->

<!-- badges: end -->

Single molecule tracking has evolved as a novel new approach
complementing genomic sequencing, it reports live biophysical properties
of molecules being investigated besides properties relating their coding
sequence.

Here we provided “sojourner” package, to address statistical and
bioinformatic needs related to the analysis and comprehension of high
throughput single molecule tracking data.

## Installation

You can install the released version from [GitHub](https://github.com/)
with:

``` r
# install.packages("devtools")
devtools::install_github("sheng-liu/sojourner", build_vignettes = TRUE)
```

## Dependencies

The packages that sojourner uses in case not installed automatically.

``` r
# plotting
install.packages('ggplot2')
install.packages('reshape2')
install.packages('dplyr')
install.packages('scales')
install.packages('gridExtra')

# curve fitting and clustering
install.packages('mixtools')
install.packages('fitdistrplus')
install.packages('nls2')
install.packages('minpack.lm')
install.packages('truncnorm')
install.packages('mclust')

# image processing - CRAN
install.packages('rtiff')
install.packages('jpeg')
install.packages('tiff')
install.packages('png')
install.packages('pixmap')

# image processing- Bioconductor
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install()
BiocManager::install(c("EBImage"))

# interfacing
install.packages('R.matlab')
install.packages('Rcpp')
install.packages("reticulate")

# GUI
install.packages('shiny')
install.packages('shinyjs')

# others
install.packages('mltools')
install.packages('sampSurf')
install.packages('sp')
install.packages('rlang')
install.packages('SDMTools')
install.packages('rowr')
```
