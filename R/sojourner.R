# sojourner help (environment)


## ----------------------------------------------------------------------------
## sojourner Roxygen help
##' @name sojourner
##' @aliases sojourner
##' @title sojourner - statistical analysis of single molecule trajectories
##' @rdname sojourner
##' @docType package
##' @description Single molecule tracking reports live biophysical properties of
##'   molecules being investigated besides their coding sequence. It has evolved
##'   as a novel approach complementing genomic sequencing. Here we provided
##'   “sojourner” package, to address statistical and bioinformatic needs
##'   related to the analysis and comprehension of high throughput single
##'   molecule tracking data.

##  @usage
## library(sojourner)

##' @details sojourner package provides statistical analysis of biophysical
##'   properties of single molecules. Current version primarily uses mean square
##'   displacement (MSD) based analysis, more functions using displacement based
##'   and hidden Markov model based method will be added in future release.
##'
##'   Current version includes: dwellTime (duration of the tracks), squareDisp
##'   (squared displacement), MSD (mean square displacement), Dcoef(diffusion
##'   coefficient), CDF(cumulative distribution function), etc., and various
##'   plotting functions for visulization. It also includes a graphical user
##'   interface (GUI) sojournerGUI(). The GUI version has the most frequently
##'   used functions but not all command line ones.
##'   
##' 
## @examples
## # not run
## # library(sojourner)
## # sojournerGUI()
##' 
##' 
##' @import ggplot2
##' @import utils
##' @import graphics
##' @import grDevices
##' @import grid
##' @import reshape2
##  @importFrom reshape2 melt

##' @importFrom stats approx coef coefficients complete.cases sigma dnorm
##' density deviance lm median nls nls.control qt sd setNames var
##' @importFrom scales cbreaks
##' @importFrom mixtools normalmixEM
##' @importFrom mixtools boot.se
##' @importFrom mixtools lambda
##' @importFrom fitdistrplus fitdist
##' @importFrom fitdistrplus denscomp
##' @importFrom nls2 nls2
##' @importFrom minpack.lm nlsLM
##' @importFrom truncnorm rtruncnorm
##' @importFrom mclust mclustBootstrapLRT
##' @importFrom gridExtra grid.arrange
##' @importFrom gridExtra marrangeGrob
##' @importFrom rtiff readTiff
##' @importFrom EBImage readImage
##' @importFrom compiler cmpfun
##' @importFrom compiler enableJIT
##' @importFrom plyr rbind.fill
##' @importFrom rlang .data
##' @importFrom shinyjs hidden
##' @importFrom shinyjs enable
##' @importFrom shinyjs disable
##' @importFrom shinyjs show
##' @importFrom shinyjs hide
##' @importFrom shinyjs delay
##' @importFrom shinyjs html
##' @importFrom shinyjs useShinyjs

## @import dplyr
## dplyr has masked intersect, setdiff, setequal, union from base and
## other packages, try to use importFrom instead of import package

##' @importFrom dplyr bind_rows
##' @importFrom dplyr select
##' @importFrom dplyr lag
##' @importFrom dplyr one_of
##' 
##' 
##' @importFrom MASS kde2d
##' @importFrom sp point.in.polygon

##' @import shiny

## required for roxygen2 to genrate rd files
NULL



