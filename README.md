
# <img src="vignettes/sojourner_logo.png" align="right" alt="" width="180" />

<!-- README.md is generated from README.Rmd. Please edit that file -->

# *Sojourner*: Statistical Analysis of Single Molecule Trajectories

<!-- badges: start -->

<!-- badges: end -->

Single molecule tracking has evolved as a novel new approach
complementing genomic sequencing, it reports live biophysical properties
of molecules being investigated besides properties relating their coding
sequence. Here we provided “*sojourner*” package, to address statistical
and bioinformatic needs related to the analysis and comprehension of
high throughput single molecule tracking data
(<https://bioconductor.org/packages/sojourner/>).

Maintainer:

*sojourner* developer \<sojourner.developer at outlook.com\>

Citation:

Liu S, Yoo S, Tang X, Sung Y, Wu C (2020). Sojourner: statistical
analysis of single molecule trajectories. R package version 1.3.0,
<https://github.com/sheng-liu/sojourner>. DOI:
10.18129/B9.bioc.sojourner.

## Installation

You can install *sojourner* as a Bioconductor package:

``` r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("sojourner")
```

Alternatively, you can also install *sojourner* directly from github.

``` r
# install.packages("devtools")
devtools::install_github("sheng-liu/sojourner", build_vignettes = TRUE, ref = "master")
```

## Documentation

*Sojourner* has a vignette contained in the package. You can view the
documentation by

``` r
browseVignettes("sojourner")
```

Alternatively, you can also access the the vignette on *sojourner*’s
website

<https://sheng-liu.github.io/sojourner/articles/sojourner-vignette.html>

## Dependencies

These are packages that *sojourner* uses in case they are not installed
automatically.

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

```
