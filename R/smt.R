# smt help (environment)
#
#
###############################################################################


## -----------------------------------------------------------------------------
## smt Roxygen help
##' @name smt
##' @aliases smt
##' @title smt-simple parameters for single molecule tracking analysis
##' @rdname smt
##' @docType package
##' @description simple analysis on single molecule tracking data using parameters based on mean square displacement (MSD).
## @usage
## smt()

##' @details smt provide a simple analysis on single molecule tracking data using parameters based on mean square displacement (MSD). Currently includes:
##' - duration of the tracks (dwellTime),
##'
##' - square displacement (squareDisp),
##'
##' - mean square displacement as a function of time (msd),
##'
##' - diffusion coefficient (Dcoef) and
##'
##' - emperical cumulative distribution function (eCDF) of MSD over time.

## @seealso

##' @import ggplot2

##' 
##' 
##' 
##' 
##' @import reshape2
## @import gridExtra
## @importFrom reshape2 melt
##' @importFrom scales cbreaks
##' @importFrom mixtools normalmixEM
##' @importFrom mixtools boot.se
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

## @import dplyr

## dplyr has masked  intersect, setdiff, setequal, union from base and other
## packages, try to use importFrom instead of import package

##' @importFrom dplyr bind_rows
##' @importFrom dplyr select
##' @importFrom dplyr lag
##' 
##' 
##' @importFrom MASS kde2d
##' @importFrom sp point.in.polygon

##' @import shiny
##' @export smtGUI

smt=function(){}

smtGUI=function(){
    appPath=system.file("myapp",package="smt")
    shiny::runApp(appPath)
f
}


# shinyjs

