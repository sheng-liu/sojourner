# sojourner help (environment)
#
#
###############################################################################


## -----------------------------------------------------------------------------
## sojourner Roxygen help
##' @name sojourner
##' @aliases sojourner sojournerGUI
##' @title sojourner-simple parameters for single molecule tracking analysis
##' @rdname sojourner
##' @docType package
##' @description simple analysis on single molecule tracking data using 
##' parameters based on mean square displacement (MSD).
##' @usage
##' sojourner()
##' sojournerGUI()


##' @details sojourner provide a simple analysis on single molecule tracking 
##' data using parameters based on mean square displacement (MSD). Currently 
##' includes:
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
##' @import utils
##' @import graphics
##' @import grDevices
##' @import grid
##' @importFrom stats approx coef coefficients complete.cases sigma dnorm
##' density deviance lm median nls nls.control qt sd setNames var
##' 
##' 

##' 
##' @import reshape2
## @import gridExtra
## @importFrom reshape2 melt
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

## @import dplyr

## dplyr has masked  intersect, setdiff, setequal, union from base and other
## packages, try to use importFrom instead of import package

##' @importFrom dplyr bind_rows
##' @importFrom dplyr select
##' @importFrom dplyr lag
##' @importFrom dplyr one_of
##' 
##' 
##' @importFrom MASS kde2d
##' @importFrom sp point.in.polygon

##' @import shiny
##' @import shinyjs
##' @export sojournerGUI
##' @return launch sojourner
sojourner=function(){}

sojournerGUI=function(){
    appPath=system.file("myapp",package="sojourner")
    shiny::runApp(appPath)
}


# shinyjs

