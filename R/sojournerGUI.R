# sojournerGUI (environment)


## ----------------------------------------------------------------------------
## sojournerGUI Roxygen help
##' @name sojournerGUI
##' @aliases sojournerGUI
##' @title sojournerGUI
##' @rdname sojournerGUI-methods
##' @docType methods
##' @description Graphical user interface (GUI) for sojourner.

##' @usage
##' sojournerGUI()

##' @details 
##' The GUI version has the most frequently used functions but not all command
##' line ones. Current version includes: MSD (mean square displacement),
##' Dcoef(diffusion coefficient), CDF(cumulative distribution function), and
##' RT(residence time) etc.
##' 
##' 
##' @examples
##' # not run
##' # library(sojourner)
##' # sojournerGUI()
##' 
## required for BiocCheck pass
##' @return A graphical user interface for sojourner.
##' 
##' @export sojournerGUI

sojournerGUI = function() {
    appPath = system.file("myapp", package = "sojourner")
    shiny::runApp(appPath)
}




